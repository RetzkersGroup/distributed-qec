Here, we present the general concepts for creating a repeated measurement floquet circuit.

# Goal

In order to reduce measurement errors, one can repeat the pairwise measurements composing a Floquet circuit.

The results can be decoded using majority vote (iner-decoder) and outer decoder (same as the normal circuit), or a full decoder that decodes everything together. The full decoder will present superior error rates and also allow for other scenarios such as 2-repetiotions (which does not make sense with majority decoder).

for the full-decoder, I want to build a floquet circuit that has the normal detectors, as well as small detectors between the measurements.


# input circuit

To be filled.

example: in/circuits/honeycomb_gidney/model=SDEM3,r=10,w=6,h=9,p=0.001,obs=H,sheared=False.stim

# Output:

## multiply the measurement

first, each time step (between two TICKs) should be precedded with r-1 new time steps. for example, r=2:

```
...
TICK
MXX 0 1 2 3
DETECTOR rec[-1]
TICK
DETECTOR rec[-8]
```

should become:

```
TICK
MXX 0 1 2 3
TICK
MXX 0 1 2 3
DETECTOR rec[-1]
TICK
DETECTOR rec[-8]
```

## add Detector

the new measurement, since it is exactly before the new one, should have the same result. so a detector is just multiplication of both. So we need to add:

```
TICK
MXX 0 1 2 3
TICK
MXX 0 1 2 3
# new detectors:
DETECTOR rec[-1] rec[-3]
DETECTOR rec[-2] rec[-4]
DETECTOR rec[-1]
TICK
DETECTOR rec[-8]
```

## Fix records (Detectors and observables):

### Detectors

since in stim target records are defined relative to the detector, we need to update the detectors below the addition:

#### Records that are not part of the repeated measurement:

Here, we just need to update the rec numbers. Like the last detector here:

```
TICK
MXX 0 1 2 3
TICK
MXX 0 1 2 3
# new detectors:
DETECTOR rec[-1] rec[-3]
DETECTOR rec[-2] rec[-4]
DETECTOR rec[-1]
TICK
# we added two measurements, so we need to update the detector by two
DETECTOR rec[-10]
```

#### Records that are part of the repeated measurement:

We ideally want to keep the decoding graph matchable, thus, we can't leave the records as they are since this might mean that a measurement that is used by 2 detectors will now be used by three, which is not matchable.\

Insteasd, we need to find which two detectors are using the original measurement as record, and point one of the to the earliest of the repeated measurement and another to the last. I am not sure yet which of the detectors should use the first repeated measurement and which should use the last. offer me options.

Also, this feature should be optional - if the original code is not matchable so there is no need for it.

### OBSERVABLE_INCLUDE

those are also relative. I am not sure yet how to fix those.