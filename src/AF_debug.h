#ifndef AF_DEBUG_H_
#define AF_DEBUG_H_
/*********************************************************************
 * DEVELOPER:
 *  - Jonathan Kim (jkm@illinois.edu)
 *
 * DESCRIPTION:
 *  Debug related handling options and functions.
 *
 */

#include <iostream>

/* =======================
 * DEBUG Options
 */
// set to 1 to show extra info for debugging for the AF tool code
#define DEBUG_TOOL  0
// set to 1 to show extra info for debugging for the AF tool code
#define DEBUG_TOOL_PARSER  0
// set to 1 to show elapsed time or performance debugging info for the AF tool code
#define DEBUG_ELAPSE_TIME 0
// set to 1 to show extra info for debugging for the io code
#define DEBUG_IO  0


/* =======================
 * Measure time functions
 */
void StartElapseTime();
void StopElapseTimeAndShow(std::string msg);


#endif // AF_DEBUG_H_
