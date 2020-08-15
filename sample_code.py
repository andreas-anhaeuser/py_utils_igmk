#!/usr/bin/python

from chronmeter_new import Chronometer

N = 3
M = 4
K = 5
J = 6
I = 7
code = None

code        # <-- [0] 
chrono = Chronometer
for n in chrono.range(N):
    code    # <-- [0, 0]

    for m in chrono.range(M):
        code    # <-- [0, 0, 0]

    code        # <-- [0, 1]
    for m in chrono.range(M):
        code    # <-- [0, 1, 0]

        for k in chrono.range(K):
            code    # <-- [0, 1, 0, 0]

        code    # <-- [0, 1, 1]

        for j in chrono.range(J):
            code    # <-- [0, 1, 1, 0]

        code    # <-- [0, 1, 2]

    code    # <-- [0, 2]

    for i in range(I):
        code    # <-- [0, 2, 0]

    code    # <-- [0, 3]

code    # <-- [1]
