/* Example: compute Sun rise, transit, set times for a given location.
 * Compile: gcc -o riseset riseset.c -I./include -L./target/release -llibre_ephemeris -lm
 */
#include "le_ephemeris.h"
#include <stdio.h>

int main() {
    // Zurich: 47.4°N, 8.5°E
    double lat = 47.4;
    double lon = 8.5;
    double jd_ut = oe_julday(2026, 6, 22, 1);

    char serr[256];
    double xx[24];

    le_set_topo(lon, lat, 0.0);

    printf("Location: %.1f°N, %.1f°E\n", lat, lon);
    printf("Date: June 22, 2026 (JD %.1f)\n\n", jd_ut);

    // RISE = 0, SET = 1, TRANSIT = 2, TRANSIT_LOWER = 3
    const char *events[] = {"Rise", "Set", "Transit", "Lower transit"};

    int bodies[] = {LE_SUN, LE_MOON};
    const char *body_names[] = {"Sun", "Moon"};

    for (int b = 0; b < 2; b++) {
        printf("--- %s ---\n", body_names[b]);
        for (int e = 0; e < 3; e++) {
            double t = le_rise_trans(jd_ut, bodies[b], 0, e, lon, lat, xx, serr);
            if (t > 0) {
                double hours = (t - jd_ut) * 24.0;
                int h = (int)hours;
                int m = (int)((hours - h) * 60);
                printf("  %s: %02d:%02d UT\n", events[e], h, m);
            } else {
                printf("  %s: below horizon\n", events[e]);
            }
        }
        printf("\n");
    }

    return 0;
}
