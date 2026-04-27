# =============================================
# RK4 desde cero — du/dt = 1.01*u, u(0) = 0.5
# =============================================

# Función f(t, u)
f(t, u) = 1.01 * u

# Condición inicial y rango de tiempo
u0     = 0.5
t0, tf = 0.0, 1.0
dt     = 0.1          # tamaño del paso

# Solución analítica para comparar
analitica(t) = 0.5 * exp(1.01 * t)

# ---- RK4 ----
function rk4(f, u0, t0, tf, dt)
    ts = t0:dt:tf
    us = zeros(length(ts))
    us[1] = u0

    for i in 1:length(ts)-1
        t, u = ts[i], us[i]
        k1 = f(t,            u)
        k2 = f(t + dt/2,     u + dt/2 * k1)
        k3 = f(t + dt/2,     u + dt/2 * k2)
        k4 = f(t + dt,       u + dt   * k3)
        us[i+1] = u + (dt/6) * (k1 + 2k2 + 2k3 + k4)
    end
    return ts, us
end

ts, us = rk4(f, u0, t0, tf, dt)

# ---- Resultados ----
println("  t      RK4        Analítica   Error")
println("─"^45)
for (t, u) in zip(ts, us)
    err = abs(u - analitica(t))
    @printf("  %.1f    %.6f   %.6f    %.2e\n", t, u, analitica(t), err)
end