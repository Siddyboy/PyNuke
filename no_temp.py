# START
#T = 0.00364

running = True

while running is True:

    inputvalid = False

    while inputvalid is False:

        print("==========================================")
        inputdata = input("Effect (%)?")

        if inputdata is not "":
            try:
                effdata = float(inputdata)
                if effdata < -98.8 or effdata > -1.0:
                    print("Input out of valid range")
                    inputvalid = False
                else:
                    inputvalid = True
            except ValueError as e:
                print("Input probably not a number...")
        elif inputdata == "":
            raise RuntimeError("No data entered: program terminated.")

    NMAX = 100

    hot = 1.0
    cold = 0.0
    TOL = 0.0000001
    n = 0

    mid = (cold + hot) / 2.0

    while n < NMAX:

        print("h=", format(hot, "0.6f"), \
              "c=", format(cold, "0.6f"), \
              "m=", format(mid, "0.6f"), \
              "eff=", format(effect(mid), "0.3f"), \
              "n=", format(n, "d"))

        if (hot - cold) / 2 < TOL:
            break
        a = effdata
        b = effect(mid)
        if a < b:
            hot = mid
        else:
            cold = mid
        mid = (cold + hot) / 2
        n += 1

    print("T = ", format(mid * 1000, "0.2f"), " mK", sep='')
