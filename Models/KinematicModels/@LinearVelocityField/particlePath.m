function [x, y] = particlePath(self, x0, y0, t, kappa, u_0, v_0)
% particlePath Integrate the linear velocity field with an Euler-Maruyama random walk term.
%
% Declaration:
%   `[x, y] = particlePath(self, x0, y0, t, kappa, u_0, v_0)`
%
% Parameters:
%   `x0`, `y0` - Initial particle positions in meters. These inputs may be
%   vectors for multiple particles and are normalized to row vectors.
%
%   `t` - Output time vector in seconds with shape `[nTimes 1]` or
%   `[1 nTimes]`.
%
%   `kappa` - Scalar diffusivity in m^2 s^-1.
%
%   `u_0`, `v_0` - Background velocities in m s^-1 used by the analytical
%   formulas.
%
% Returns:
%   `x`, `y` - Position arrays in meters with shape `[length(t) nParticles]`.
arguments
    self (1,1) LinearVelocityField
    x0 {mustBeNumeric}
    y0 {mustBeNumeric}
    t {mustBeNumeric}
    kappa (1,1) double
    u_0 (1,1) double
    v_0 (1,1) double
end

x0 = reshape(x0, 1, []);
y0 = reshape(y0, 1, []);
t = reshape(t, [], 1);

deltaT = t(2) - t(1);
nTimes = length(t);

randAmp = sqrt(deltaT * 2 * kappa);
dX = randAmp * randn(nTimes, length(x0));
dY = randAmp * randn(nTimes, length(y0));

sigma2 = self.sigma_n * self.sigma_n + self.sigma_s * self.sigma_s;
s2 = sigma2 - self.zeta * self.zeta;
s = sqrt(s2);

x = ones(length(t), 1) * x0;
y = ones(length(t), 1) * y0;

for iTime = 2:length(t)
    xn = x(iTime - 1,:);
    yn = y(iTime - 1,:);

    if sigma2 > 0 && s2 == 0
        x(iTime,:) = (1 + self.sigma_n * deltaT / 2) * xn + ((self.sigma_s - self.zeta) * deltaT / 2) * yn;
        y(iTime,:) = ((self.sigma_s + self.zeta) * deltaT / 2) * xn + (1 - self.sigma_n * deltaT / 2) * yn;
    elseif s2 < 0
        s = sqrt(self.zeta * self.zeta - sigma2);
        cos_t = cos(s * deltaT / 2);
        sin_t = sin(s * deltaT / 2);

        x(iTime,:) = (cos_t + (self.sigma_n / s) * sin_t) * xn + ((self.sigma_s - self.zeta) / s) * sin_t * yn;
        y(iTime,:) = ((self.sigma_s + self.zeta) / s) * sin_t * xn + (cos_t - (self.sigma_n / s) * sin_t) * yn;

        x(iTime,:) = x(iTime,:) + (2 / s^2) * ((s * sin_t + self.sigma_n * (1 - cos_t)) * u_0 + (self.sigma_s - self.zeta) * (1 - cos_t) * v_0);
        y(iTime,:) = y(iTime,:) + (2 / s^2) * ((self.sigma_s + self.zeta) * (1 - cos_t) * u_0 + (s * sin_t - self.sigma_n * (1 - cos_t)) * v_0);
    elseif s2 > 0
        cosh_t = cosh(s * deltaT / 2);
        sinh_t = sinh(s * deltaT / 2);

        x(iTime,:) = (cosh_t + (self.sigma_n / s) * sinh_t) * xn + ((self.sigma_s - self.zeta) / s) * sinh_t * yn;
        y(iTime,:) = ((self.sigma_s + self.zeta) / s) * sinh_t * xn + (cosh_t - (self.sigma_n / s) * sinh_t) * yn;

        x(iTime,:) = x(iTime,:) + (2 / s^2) * ((s * sinh_t + self.sigma_n * (cosh_t - 1)) * u_0 + (self.sigma_s - self.zeta) * (cosh_t - 1) * v_0);
        y(iTime,:) = y(iTime,:) + (2 / s^2) * ((self.sigma_s + self.zeta) * (cosh_t - 1) * u_0 + (s * sinh_t - self.sigma_n * (cosh_t - 1)) * v_0);
    else
        x(iTime,:) = xn;
        y(iTime,:) = yn;
    end

    x(iTime,:) = x(iTime,:) + dX(iTime,:);
    y(iTime,:) = y(iTime,:) + dY(iTime,:);
end
