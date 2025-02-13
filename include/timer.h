/* Header file for C-compatible variables/functions in timer.fF90 */

/* Public variables within the timer module */
extern double timer_qtime;     /* Walltime of batch queue in seconds */
extern double timer_closetime; /* Time needed for clean shutdown     */ 

/* Initialise timer */
void timer_init(void);

/* Check if we are within timer_closetime of the maximum runtime */
/* returns 0 if not.                                             */
int timer_check_continuation(void);

/* Function to check how much time has elapsed since the start of the run */
double timer_elapsed_time(void);

