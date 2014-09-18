#ifndef PROFILERC_H
#define PROFILERC_H
extern "C" {
void* profilerNew(char* name);
void profilerReset(void* profiler, char* name);
void profilerStart(void* profiler, char* name);
void profilerDeclare(void* profiler, char* name);
void profilerStop(void* profiler, char* name, long operations=0);
char* profilerStr(void* profiler);
void profilerStrSubroutine(void*profiler, char* result, int maxResult);
}

#endif // PROFILERC_H
