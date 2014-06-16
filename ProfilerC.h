#ifndef PROFILERC_H
#define PROFILERC_H
extern "C" {
void* profilerNew(char* name);
void profilerReset(void* profiler, char* name);
void profilerStart(void* profiler, char* name);
void profilerStop(void* profiler, char* name);
char* profilerStr(void* profiler);
}

#endif // PROFILERC_H
