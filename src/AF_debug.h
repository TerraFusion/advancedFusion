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

//-----------------------
// DEBUG Options
// 1 will show extra info for debugging for the AF tool code
#define DEBUG_TOOL  0
// 1 will show elapsed time or performance debugging info for the AF tool code
#define DEBUG_ELAPSE_TIME 0

//-------------------------
// Measure time functions
void StartElapseTime();
void StopElapseTimeAndShow(std::string msg);


#endif // AF_DEBUG_H_
