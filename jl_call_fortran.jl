lib = ("./mod_diag_functions.so")

FT=Float64

"""
  calc_rh(p, t, qv)

# Arguments

- `p::FT`: pressure (Pa)
- `t::FT`: temperature (K)
- `qv::FT`: water vapor mixing ratio (kg/kg)

# Examples

```julia
calc_rh(1000_00, 20. + K0, 3.0/1e3)
```
"""
function calc_rh(p, t, w)
  # Constants
  pq0 = 379.90516
  a2 = 17.2693882
  a3 = 273.16
  a4 = 35.86
  rhmin = 1.0

  # algorithms adapted from WRFPOST
  q = w / (1.0 + w)
  qs = pq0 / p * exp(a2 * (t - a3) / (t - a4))
  rh = 100.0 * q / qs
  clamp(rh, rhmin, 100.0)
end


K0 = 273.15
tK = 20 + K0
# es: hPa
cal_es(tK::FT) = @ccall lib.__diag_functions_MOD_calc_es(tK::Ref{FT})::FT

cal_Tdew(tC::FT, rh::FT) = 
  @ccall lib.__diag_functions_MOD_calc_dewpoint(tC::Ref{FT}, rh::Ref{FT})::FT


# cal_Tdew(tC::FT, rh::FT) =
theta(t::FT, p::FT) = 
  @ccall lib.__diag_functions_MOD_theta(t::Ref{FT}, p::Ref{FT})::FT

"""
  theta_se(tK::FT, p::FT, rh::FT, w::FT)

- tK: K
- p: hPa
- rh: %
- w: g/kg
"""
theta_se(tK::FT, p::FT, rh::FT, w::FT) =
  @ccall lib.__diag_functions_MOD_thetae(
    tK::Ref{FT}, p::Ref{FT}, rh::Ref{FT}, w::Ref{FT})::FT

# 每个人算出来的都不太一样
theta_se(20 +K0, 1000.0, 20.0, 0.002884103) - K0

# atm = 1000.0 #* 100.0 # hPa
# theta(20 + K0, atm)  - K0

cal_Tdew(20., 80.)
cal_es(tK)
