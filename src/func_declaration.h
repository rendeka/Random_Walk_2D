#pragma once

void task_1(int trajectory[n][2], int step, bool dirOptions[4], double walkLength[N], double avgWalkLength);
void make_step(int trj[n][2], int dir, int i, int step);
bool step_dont_cross(int trj[n][2], int* i, int step, bool dirOpts[4], int* prevDir);
bool is_not_dead_end(bool dirOpts[4]);
int next_direction(bool dirOpts[4]);
double norm(int x, int y);
