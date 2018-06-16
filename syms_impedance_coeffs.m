syms Rr Xr Rs Xs Rm Xm s

Zr = Rr/s + 1i*Xr
Zs = Rs + 1i*Xs
Zm = 1./(1./Rm + 1./(1i*Xm))

Zgen = Zs + 1/(1/Zm + 1/Zr)

pretty(simplifyFraction(Zgen))