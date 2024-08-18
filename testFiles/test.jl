using Trapz: trapz
vx=0:0.1:1
vy=(x->x^2).(vx)
trapz(vx, vy) / (vx[end] - vx[1]) #
mean(vy)