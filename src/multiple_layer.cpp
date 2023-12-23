#include "multiple_layer.h"

/**
 * Simulate with a multpile-layer tissue media:
 * 
 *             n0 = 1
 *    -----------------------------
 * 
 *   t1 = 0.05mm     n1 = 1.4    mu_a1 = 37 cm-1, mu_s1 = 480 cm-1, g1 = 0.79
 * 
 *    -----------------------------
 * 
 *   t2 = 2mm        n2 = 1.5    mu_a2 = 2.2 cm-1, mu_s2 = 220 cm-1, g2 = 0.79
 * 
 *    -----------------------------
 * 
 *             n0 = 1
 * 
 * Parameters:
 * 
 *    1. delta r = delta z = 0.025 mm
 *    2. variable weight photons
 *    3. Henyey-Greenstein phase function
 * 
*/