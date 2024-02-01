#ifndef CONFIG_H
#define CONFIG_H

//#define SKIP_OPTIMAL
//#define SKIP_NAIVE

// Measurement testing params
#define MIN_DIM                  (3)
#define STEP                     (100)
#define MAX_DIM                  (1003)

// Axes params
#define Y_MAX_VAL                (1)
#define X_MAX_VAL                (MAX_DIM)

// Velicina kanvasa
#define CANVAS_WIDTH             (MAX_DIM)
#define CANVAS_HEIGHT            (MAX_DIM)

#define BROJ_SLUCAJNIH_OBJEKATA  (20)
#define DUZINA_PAUZE             (500)

/* Ne prevelika tolerancija na numericku gresku */
#define EPS                      (1e-6)
#define EPSf                     (1e-6f)

#endif // CONFIG_H
