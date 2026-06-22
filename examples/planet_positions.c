/* Example: compute planetary positions for today.
 * Compile: gcc -o planet_positions planet_positions.c -I./include -L./target/release -llibre_ephemeris -lm
 */
#include "le_ephemeris.h"
#include <stdio.h>
#include <math.h>

int main() {
    double jd_ut = le_julday(2026, 6, 22.5, 1);
    double xx[24];
    char serr[256];

    int flags = LE_FLG_XYZ | LE_FLG_J2000
        | LE_FLG_NOABERR | LE_FLG_NOGDEFL | LE_FLG_NOBIRR | LE_FLG_NONUT;

    printf("JD = %.1f (June 22, 2026)\n\n", jd_ut);

    int planets[] = {LE_MERCURY, LE_VENUS, LE_MARS, LE_JUPITER, LE_SATURN,
                     LE_URANUS, LE_NEPTUNE, LE_PLUTO};
    const char *names[] = {"Mercury", "Venus", "Mars", "Jupiter", "Saturn",
                          "Uranus", "Neptune", "Pluto"};

    for (int i = 0; i < 8; i++) {
        int rc = le_calc_ut(jd_ut, planets[i], flags, xx, serr);
        if (rc != LE_OK) {
            printf("%-8s: error %d\n", names[i], rc);
            continue;
        }
        double dist = sqrt(xx[LE_X]*xx[LE_X] + xx[LE_Y]*xx[LE_Y] + xx[LE_Z]*xx[LE_Z]);
        printf("%-8s: (%.4f, %.4f, %.4f)  dist=%.4f AU\n",
               names[i], xx[LE_X], xx[LE_Y], xx[LE_Z], dist);
    }

    // Moon position (geocentric)
    int rc = le_calc_ut(jd_ut, LE_MOON, flags, xx, serr);
    if (rc == LE_OK) {
        double dist = sqrt(xx[LE_X]*xx[LE_X] + xx[LE_Y]*xx[LE_Y] + xx[LE_Z]*xx[LE_Z]);
        printf("\nMoon (geocentric): (%.4f, %.4f, %.4f)  dist=%.4f AU = %.1f km\n",
               xx[LE_X], xx[LE_Y], xx[LE_Z], dist, dist * 149597870.7);
    }

    return 0;
}
