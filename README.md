# PQRST Monitor Firmware

Single-lead ECG monitor for ESP32. An AD8232 front end is sampled at 250 Hz and the signal is filtered. Beats are detected with a Pan-Tompkins QRS detector, each beat is measured, and the results go out over UART together with the filtered waveform.

```
dsp/              portable C99 signal processing (no malloc, no platform headers)
  ecg_filter      HPF 0.5 Hz, notch 60 Hz, LPF 40 Hz
  qrs_detector    Pan-Tompkins with adaptive thresholds, refractory, searchback
  ecg_delineate   QRS on/off, P onset, T end, baseline, R and ST amplitude
  ecg_monitor     per-sample pipeline: filter -> detect -> delineate -> beat
  ecg_report      UART line protocol
firmware/esp32/   ESP-IDF 5.x project (timer-driven sampling, tasks, UART)
test/             host test harness on synthetic ECG with known fiducials
```

## Wiring (ESP32 DevKitC + AD8232 breakout)

| AD8232 | ESP32 | Notes |
|---|---|---|
| 3.3V | 3V3 | do not power from 5 V |
| GND | GND | |
| OUTPUT | GPIO36 (ADC1_CH0, VP) | must be an ADC1 pin |
| LO+ | GPIO34 | input-only pin, AD8232 drives it |
| LO- | GPIO35 | |
| SDN | 3V3 | or a GPIO to power down the front end |

Electrodes: RA (red), LA (yellow), RL/right-leg drive (green), Lead I or II placement.

All pins, the front-end gain (1100) and the UART settings live in `firmware/esp32/main/board.h`. Set `ECG_LO_P_GPIO` / `ECG_LO_N_GPIO` to `-1` if the lead-off pins are not connected.

The AD8232 output sits at mid-rail, 1.65 V. With a gain of 1100, the input range before clipping is about +0.7 mV on the original ESP32, whose calibrated 12 dB range is linear to ~2.45 V. On chips with curve-fitting calibration it is about ±1.3 mV. Reduce the gain on the breakout if R waves are taller than that.

## Build

Firmware (ESP-IDF 5.1 or newer):

```
cd firmware/esp32
idf.py set-target esp32
idf.py build flash monitor
```

Host tests (any C99 compiler):

```
make -C test run
./test/test_ecg dump.csv         # also writes raw/filtered traces
```

## Signal chain

```
ADC -> mV --> HPF 0.5 Hz --> notch 60 Hz --> LPF 40 Hz --+--> UART stream
              linear phase   Q 30, 2 Hz BW   Butterworth  |
                                             4th order    +--> Pan-Tompkins --> delineation --> beat
```

- **High-pass 0.5 Hz.** The default is a linear-phase design: the input delayed by L−1 samples, minus two cascaded length-L moving averages, with L = 286. That gives −3 dB at 0.5 Hz, under 0.42 dB of passband ripple, and an O(1) cost per sample. An IIR high-pass at this corner has a non-linear phase response, which drags the ST segment down by about 0.05 mV and shortens QT (see Verification). The linear-phase filter adds a constant 1.14 s delay. Build with `-DECG_HPF_LINEAR_PHASE=0` for a 2nd-order Butterworth IIR with negligible delay, if display latency matters more than ST accuracy.
- **Notch 60 Hz.** Biquad notch, more than 60 dB rejection at 60 Hz. Set `ECG_NOTCH_HZ` to `50.0f` for 50 Hz mains.
- **Low-pass 40 Hz.** 4th-order Butterworth (two biquads).
- Filter state is primed from the first sample, so the 1.65 V front-end offset doesn't cause a start-up transient.

All corners and detector/delineation parameters are in `dsp/ecg_config.h`, and each one can be overridden with `-D`.

## QRS detection

The Pan-Tompkins pipeline runs on the filtered signal:
- 5–15 Hz band-pass
- five-point derivative
- squaring
- 150 ms moving-window integration

Decision rules:
- **Training.** The first 2 s set the initial signal and noise peak levels (SPKI/NPKI).
- **Thresholds.** THR1 = NPKI + 0.25 (SPKI − NPKI) and THR2 = THR1 / 2, updated after every MWI peak.
- **Refractory period.** 200 ms.
- **T-wave check.** A peak within 360 ms of the previous QRS whose maximum slope is under half that QRS's slope is classified as a T wave.
- **Searchback.** If no beat is found for 166% of the mean RR (last 8 beats), the largest peak above THR2 is accepted.
- **Amplitude recovery.** SPKI is halved after 2.5 s without a beat, so the detector recovers after an amplitude drop.

## Per-beat measurements

| Field | Method |
|---|---|
| R index | dominant QRS deflection in the filtered signal; polarity is tracked so RS complexes don't flip |
| RR, HR | consecutive R peaks, refined to sub-sample by parabolic interpolation; HR = 60000 / RR |
| QRS duration | onset/offset where the mean \|slope\| over 20 ms falls below max(4% of peak QRS slope, EMG noise floor) |
| Baseline | mean of the flattest 20 ms in the PR segment |
| PR interval | P onset (tangent at the steepest P upstroke meets baseline) to QRS onset |
| QT interval | QRS onset to T end (tangent at the steepest T downstroke meets baseline) |
| R amplitude | max of QRS minus baseline, mV |
| ST level | J + 60 ms minus baseline, mV (J + 40 ms when HR > 120 bpm) |

A beat is reported once its T wave has been seen. That happens when the next QRS is detected, or after 0.7 RR + 0.3 s, whichever comes first. Measurements that fail sanity limits are reported as empty fields rather than guessed:
- P wave or T wave not found
- PR outside 60–400 ms
- QT outside 150–700 ms
- QRS outside 30–220 ms

## UART protocol

The protocol is ASCII at 115200 8N1, one record per line ending in `\r\n`. `CS` is the XOR of the bytes between `$` and `*`, as two hex digits. Lines that fail the checksum should be discarded; ROM boot messages on the same port are filtered out this way.

```
$I,250,1,285*5C                                    fs Hz, protocol version, filter delay (samples)
$S,4321,-1234*7E                                   sample index, filtered ECG in uV
$B,123456,72.3,830,92,156,384,1.234,-0.050*4C      beat
$E,99,LOFF*46                                      event
```

Beat fields: `n, hr_bpm, rr_ms, qrs_ms, pr_ms, qt_ms, r_mV, st_mV`. Here `n` is the R-peak index in the same timeline as `$S`. Missing values are empty, e.g. `$B,7,,,88,,,0.980,*5A`.

| Event | Meaning |
|---|---|
| `LOFF` / `LON` | leads off / reconnected (debounced 100 ms; reconnect restarts filters and training) |
| `OVR,<k>` | sampler woke late k times; the samples were replicated to keep the 250 Hz clock |
| `DROP,<k>` | k samples lost because the DSP task fell behind (short gaps are bridged) |
| `TXD,<k>` | k sample lines skipped because the UART link was saturated (beat lines are never dropped) |

The sample stream uses about 5.6 kB/s, roughly half of the 115200 baud link.

Minimal host-side parser:

```python
import functools

def parse(line):
    body, _, cs = line.strip().lstrip("$").partition("*")
    if not cs or int(cs, 16) != functools.reduce(lambda a, c: a ^ ord(c), body, 0):
        return None
    tag, *fields = body.split(",")
    return tag, [float(f) if f else None for f in fields]
```

## Firmware timing

- A GPTimer alarm at 250 Hz wakes a high-priority task on core 1 through a task notification. That task does the ADC read (oneshot driver, which is not ISR-safe) and the lead-off read, then posts the sample to a 64-entry queue without waiting. It contains no delays and no UART or DSP work.
- The DSP task on core 0 drains the queue, runs the pipeline and writes to the UART TX ring buffer. Sample lines are dropped rather than blocking when the buffer is full.
- ADC readings are converted with eFuse calibration (curve or line fitting, depending on the chip) and divided by the front-end gain, giving electrode-referred mV.

## Verification

`test/test_ecg` builds lead-II beats from Gaussian P, Q, R, S and T waves with known fiducial points. T timing follows Bazett with QTc 400 ms. Except for the clean case, the signal is corrupted with:
- 0.25 Hz and 0.07 Hz baseline wander (0.3–0.5 mV)
- 0.3 mV of 60 Hz mains
- white EMG noise
- a 1.65 V front-end offset
- 12-bit ADC quantisation

Results with the default configuration (mean absolute error per beat):

| Scenario | Se / PPV | QRS | PR | QT | ST |
|---|---|---|---|---|---|
| clean, 72 bpm | 100 / 100 % | 7.3 ms | 8.9 ms | 8.3 ms | 0.002 mV |
| noisy, 72 bpm | 100 / 100 % | 7.3 ms | 8.5 ms | 13.8 ms | 0.024 mV |
| brady, 45 bpm | 100 / 100 % | 7.5 ms | 7.9 ms | 16.7 ms | 0.024 mV |
| tachy, 150 bpm | 100 / 100 % | 7.7 ms | 15.4 ms * | 16.4 ms | 0.023 mV |
| EMG 0.05 mV rms | 100 / 100 % | 11.8 ms | 20.1 ms | 24.9 ms | 0.031 mV |
| ST +0.20 mV | 100 / 100 % | 8.6 ms | 7.7 ms | 17.2 ms | 0.018 mV (reads +0.196) |
| ST −0.15 mV | 100 / 100 % | 10.2 ms | 7.1 ms | 11.6 ms | 0.016 mV (reads −0.146) |
| irregular RR ±150 ms | 100 / 100 % | 8.9 ms | 10.0 ms | 18.8 ms | 0.030 mV |

\* At 150 bpm the P wave merges with the preceding T wave. P is found in 46% of beats and is reported empty otherwise.

With the IIR high-pass, detection is unchanged but ST reads about −0.05 mV on a flat ST segment, and a +0.20 mV elevation reads +0.115 mV.

The filter test checks −3 dB at 0.5 Hz and 40 Hz, more than 60 dB rejection at 60 Hz, ≤0.5 dB ripple from 1 to 30 Hz, and no start-up transient from a DC offset. The protocol test checks formatting, empty fields, checksums and buffer overflow.

## Limitations

- Validation so far uses synthetic signals only. Before relying on the numbers, run the pipeline against annotated recordings: MIT-BIH Arrhythmia for detection and the QT Database for delineation, resampled to 250 Hz. Then bench-test with an ECG simulator through the actual AD8232.
- The design is single-lead. P/T morphology varies by lead, and inverted or biphasic T waves are handled through the |deviation| search but have not been specifically validated.
- Wide or ectopic beats (PVCs) are detected and measured but not classified.
- With the default high-pass, the stream and beat reports lag the input by 1.14 s. Beat reports arrive a further ~0.6–0.9 s after their R peak.
- This is not a certified medical device.
