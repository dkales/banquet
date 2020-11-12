#!/usr/bin/env python3
import sys

SCALING_FACTOR = 3600.0 * 1000
SCALING_FACTOR = 1000

with open(sys.argv[1], "r") as f:
    content = f.read()

testruns = content.split("Instance: ")
if len(testruns) > 1:
    testruns.pop(0)

for test in testruns:
    lines = test.splitlines()
    # first line is instance:
    print(lines[0])
    lines.pop(0)
    # second line is header:
    print(lines[0])
    lines.pop(0)

    count = 0
    keygen, sign, ver, size, ser, deser = 0, 0, 0, 0, 0, 0

    for line in lines:
        if len(line.strip()) == 0:
            continue
        vals = line.strip().split(",")
        keygen += int(vals[0])
        sign += int(vals[1])
        ver += int(vals[2])
        size += int(vals[3])
        ser += int(vals[4])
        deser += int(vals[5])
        count += 1

    keygen = (keygen / SCALING_FACTOR) / count
    sign = (sign / SCALING_FACTOR) / count
    ver = (ver / SCALING_FACTOR) / count
    size = float(size) / 1024 / count
    ser = (ser / SCALING_FACTOR) / count
    deser = (deser / SCALING_FACTOR) / count
    print("{:.2f},{:.2f},{:.2f},{:.3f},{:.2f},{:.2f}".format(
        keygen, sign, ver, size, ser, deser))
    print("-"*80)
