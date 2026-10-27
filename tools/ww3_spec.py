"""Reader for WAVEWATCH III ASCII point spectra (the gfswave.<station>.spec files on NOMADS)."""
import numpy as np


def read_spec(path):
    """Return (freqs Hz [nk], dirs rad as written [nth], times [str], E [nt, nth, nk] m^2/Hz/rad, meta [dict per time])."""
    with open(path) as f:
        head = f.readline().split("'")
        nk, nth, nsta = (int(x) for x in head[2].split())
        tok = []
        def take(n):
            while len(tok) < n:
                tok.extend(f.readline().split())
            out = [float(x) for x in tok[:n]]
            del tok[:n]
            return out
        freqs = np.array(take(nk))
        dirs = np.array(take(nth))
        times, specs, meta = [], [], []
        while True:
            line = f.readline()
            if not line.strip():
                break
            t = line.split()
            times.append(f"{t[0][:4]}-{t[0][4:6]}-{t[0][6:8]}T{t[1][:2]}:00")
            for _ in range(nsta):
                st = f.readline()
                q = st.split("'")
                vals = q[2].split()
                meta.append(dict(name=q[1].strip(), lat=float(vals[0]), lon=float(vals[1]), depth=float(vals[2]),
                                 u10=float(vals[3]), udir=float(vals[4])))
                # frequency varies fastest: ((E(ik,ith), ik=1..nk), ith=1..nth)
                specs.append(np.array(take(nk * nth)).reshape(nth, nk))
    return freqs, dirs, times, np.array(specs), meta


def band_widths(freqs):
    """WW3 frequency bin widths for a geometric frequency grid."""
    x = freqs[1] / freqs[0]
    return freqs * (x - 1 / x) / 2


if __name__ == "__main__":
    import sys
    f, d, t, E, m = read_spec(sys.argv[1])
    df, dth = band_widths(f), 2 * np.pi / len(d)
    hs = 4 * np.sqrt((E * df[None, None, :]).sum(axis=(1, 2)) * dth)
    print(len(f), "freqs", len(d), "dirs", len(t), "times", t[0], "->", t[-1], m[0])
    for i in range(0, min(len(t), 13), 3):
        print(t[i], f"Hs {hs[i]:.2f} m")
