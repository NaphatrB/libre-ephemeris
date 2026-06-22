// MIT License
//
// Copyright (c) 2026 Libre Ephemeris Contributors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


/// Evaluate Chebyshev polynomial series and its derivative.
///
/// `coeffs`: array of N Chebyshev coefficients
/// `t`: normalized time in [-1, 1]
/// `n`: number of coefficients to use
///
/// Returns (position, velocity_derivative).
/// The velocity needs to be scaled by 2/segment_length to get actual velocity.
pub fn chebyshev_eval(coeffs: &[f64], t: f64, n: usize) -> (f64, f64) {
    // Clenshaw's algorithm for Chebyshev series:
    // T_k(t) = cos(k * acos(t))
    // The recurrence: T_0 = 1, T_1 = t, T_k = 2t * T_{k-1} - T_{k-2}
    // Position = sum_{k=0}^{n-1} coeffs[k] * T_k(t)
    // Velocity derivative (w.r.t. t) = sum_{k=0}^{n-1} coeffs[k] * T'_k(t)
    // where T'_k satisfies: T'_0 = 0, T'_1 = 1, T'_k = 2t * T'_{k-1} + 2 * T_{k-1} - T'_{k-2}

    if n == 0 {
        return (0.0, 0.0);
    }

    let mut b_k1 = 0.0;
    let mut b_k = 0.0;

    for k in (0..n).rev() {
        let b_k2 = b_k1;
        b_k1 = b_k;
        b_k = coeffs[k] + 2.0 * t * b_k1 - b_k2;
    }

    let pos = b_k - t * b_k1;

    let vel = chebyshev_derivative(coeffs, t, n);

    (pos, vel)
}

fn chebyshev_derivative(coeffs: &[f64], t: f64, n: usize) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    // T'_k(t) = k * U_{k-1}(t)
    // U_0 = 1, U_1 = 2t, U_k = 2t * U_{k-1} - U_{k-2}
    // vel = Σ coeffs[k] * k * U_{k-1}(t)
    let mut u_prev = 1.0;  // U_0
    let mut u_curr = 2.0 * t; // U_1
    let mut vel = coeffs[1];   // k=1: T'_1 = 1 = 1 * U_0
    for k in 2..n {
        // u_curr is U_{k-1}, vel should use U_{k-1} for T'_k
        vel += coeffs[k] * k as f64 * u_curr;
        let u_next = 2.0 * t * u_curr - u_prev;
        u_prev = u_curr;
        u_curr = u_next;
    }
    vel
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chebyshev_constant() {
        // If only T_0 coefficient is non-zero, result should be constant
        let coeffs = [5.0, 0.0, 0.0];
        let (pos, vel) = chebyshev_eval(&coeffs, 0.5, 3);
        assert!((pos - 5.0).abs() < 1e-12);
        assert!((vel).abs() < 1e-12);
    }

    #[test]
    fn test_chebyshev_linear() {
        // T_1(t) = t
        let coeffs = [0.0, 1.0, 0.0];
        let (pos, vel) = chebyshev_eval(&coeffs, 0.5, 3);
        assert!((pos - 0.5).abs() < 1e-12);
        // derivative of T_1 is 1
        assert!((vel - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_chebyshev_quadratic() {
        // T_2(t) = 2t^2 - 1
        let coeffs = [0.0, 0.0, 1.0];
        let (pos, vel) = chebyshev_eval(&coeffs, 0.5, 3);
        assert!((pos - (2.0 * 0.25 - 1.0)).abs() < 1e-12);
        // derivative of T_2 w.r.t. t is 4t; at t=0.5, vel = 2.0
        // Note: the derivative returned is d(pos)/d(t) where t ∈ [-1, 1]
        assert!((vel - 2.0).abs() < 1e-6);
    }
}
