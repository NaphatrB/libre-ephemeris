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

/// Heliacal phenomena: first/last visibility of celestial bodies.
///
/// Computes the dates when a body becomes visible or invisible
/// due to its proximity to the Sun. Based on standard criteria
/// from Schaefer (1991) and Meeus (1998).

use crate::constants;

/// Heliacal event type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HeliacalEvent {
    FirstMorningRise,
    LastEveningSet,
    FirstEveningRise,
    LastMorningSet,
}

/// Compute the approximate arc of vision for heliacal visibility.
///
/// The arc of vision is the angular distance from the Sun to the body
/// along the horizon, measured at the moment when the body is at a given
/// altitude. A body is visible when the arc of vision exceeds a threshold
/// that depends on the body's magnitude and atmospheric conditions.
///
/// Arguments:
/// - `elongation`: angular separation from Sun (radians)
/// - `altitude`: body altitude (radians)
/// - `sun_altitude`: Sun altitude (radians, negative = below horizon)
///
/// Returns: arc of vision in radians, or None if not applicable.
pub fn arc_of_vision(elongation: f64, altitude: f64, sun_altitude: f64) -> Option<f64> {
    if sun_altitude >= 0.0 {
        return None;
    }
    let cos_aov = (sun_altitude.sin() * altitude.sin()
        + sun_altitude.cos() * altitude.cos() * elongation.cos())
    .clamp(-1.0, 1.0);
    Some(cos_aov.acos())
}

/// Check if a body is visible given its magnitude and arc of vision.
///
/// Uses the Schaefer (1991) visibility criterion:
/// A body is visible if the arc of vision exceeds a threshold
/// that depends on its visual magnitude.
pub fn is_visible(magnitude: f64, aov_rad: f64) -> bool {
    let aov_deg = aov_rad * constants::LE_RAD;
    let threshold = match magnitude {
        m if m < -1.0 => 5.0,
        m if m < 0.0 => 8.0,
        m if m < 1.0 => 12.0,
        m if m < 2.0 => 15.0,
        m if m < 3.0 => 20.0,
        m if m < 4.0 => 25.0,
        m if m < 5.0 => 30.0,
        _ => 35.0,
    };
    aov_deg > threshold
}

/// Estimate the date of a heliacal event by searching around a given date.
///
/// This is a simplified search that checks visibility at sunrise/sunset
/// for a range of dates. A full implementation would use iterative
/// refinement with accurate ephemeris calculations.
///
/// Arguments:
/// - `jd_approx`: approximate Julian day to start search
/// - `body`: planet index
/// - `event`: type of heliacal event
/// - `latitude`: observer latitude (degrees)
/// - `longitude`: observer longitude (degrees)
///
/// Returns: estimated Julian day of the event, or None if not found.
pub fn estimate_heliacal_event(
    jd_approx: f64,
    _body: i32,
    _event: HeliacalEvent,
    _latitude: f64,
    _longitude: f64,
) -> Option<f64> {
    // Placeholder: returns the approximate date.
    // Full implementation requires computing body and Sun positions
    // at sunrise/sunset for each candidate date and checking visibility.
    Some(jd_approx)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arc_of_vision_sun_above_horizon() {
        let aov = arc_of_vision(0.5, 0.3, 0.1);
        assert!(aov.is_none(), "Sun above horizon should return None");
    }

    #[test]
    fn test_arc_of_vision_sun_below_horizon() {
        let aov = arc_of_vision(0.5, 0.3, -0.1);
        assert!(aov.is_some(), "Sun below horizon should return Some");
        assert!(aov.unwrap() > 0.0);
    }

    #[test]
    fn test_is_visible_bright_object() {
        assert!(is_visible(-2.0, 0.2));
    }

    #[test]
    fn test_is_visible_dim_object() {
        assert!(!is_visible(6.0, 0.2));
    }
}
