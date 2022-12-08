#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <omp.h>
#include "settings.h"
#include "random_generators.h"
#include "func_declaration.h"

int main() 
{
	int trajectory[n][2] = { 0 };
	int step = 1;
	bool dirOptions[4] = { true, true, true, true }; // each element of this array represents a sigle direction, if the value is true then is the corresponding direction allowed
	double walkLength[N] = { 0 }; // here we store the length of each walk
	double avgWalkLength = 0; // here we store the final Euclidean distance from the origin

	if (task == 1) // simple walk
	{
	// /* // simple walk for multiple walk lengths
		int nn;
		double nnResults[powers] = { 0 };
		int nnIdx = 0;

		for (int j = 0; j < powers; j++) // we probably should create subsets of trajectory with the length of: (nn += nIncrease) but we would need to use dynamic memory for that.. maybe later
		{
			nn = std::pow(2, j + minPowerOfTwo);

			omp_set_num_threads(N_THREAD);
			#pragma omp parallel for
			for (int iPar = 0; iPar < N; iPar++)
			{
				for (int i = 0; i < nn - 1; i++)
				{
					int dir = next_direction(dirOptions);
					make_step(trajectory, dir, i, step);
				}
				walkLength[iPar] = norm(trajectory[nn - 1][0], trajectory[nn - 1][1]);
			}
			for (int i = 0; i < N; i++)
			{
				avgWalkLength += walkLength[i];
				walkLength[i] = 0; // reset for next iteration just to be sure
			}
			avgWalkLength = avgWalkLength / N;

			for (int i = 0; i < nn; i++)
			{
				trajectory[i][0] = 0;
				trajectory[i][1] = 0;
			}
			nnResults[nnIdx] = avgWalkLength;
			nnIdx++;
		}
		std::ofstream myfile;
		std::string fileName = "output/task_1/log_length_of_trj_evolution.dat";
		myfile.open(fileName);
		
		for (int j = 0; j < powers; j++)
			myfile << std::to_string(std::log(2) * ((long long int)j + minPowerOfTwo)) + "\t" + std::to_string(std::log(nnResults[j])) + "\n";

		myfile.close();
	// */

	/* // single simple walk
		#pragma omp parallel for
		for (int iPar = 0; iPar < N; iPar++)
		{
			for (int i = 0; i < n - 1; i++)
			{
				int dir = next_direction(dirOptions);
				make_step(trajectory, dir, i, step);
			}
			walkLength[iPar] = norm(trajectory[n - 1][0], trajectory[n - 1][1]);
		}
		for (int i = 0; i < N; i++)
			avgWalkLength += walkLength[i];
		avgWalkLength = avgWalkLength / N;
		std::cout << avgWalkLength << std::endl;
	*/
		std::cout << "Task 1 completed" << std::endl;
	}

	if (task == 2) // not returning walk
	{
		int dir;
		int prevDir = -1;
		int nn;
		double nnResults[powers] = { 0 };
		int nnIdx = 0;

		for (int j = 0; j < powers; j++)
		{
			nn = std::pow(2, j + minPowerOfTwo);

			omp_set_num_threads(N_THREAD);
			#pragma omp parallel for
			for (int iPar = 0; iPar < N; iPar++)
			{
				for (int i = 0; i < nn - 1; i++)
				{
					if (prevDir != -1)
						dirOptions[prevDir] = false;

					dir = next_direction(dirOptions);
					make_step(trajectory, dir, i, step);

					for (int q = 0; q < 4; q++)
						dirOptions[q] = true;

					prevDir = (dir + 2) % 4;
					walkLength[iPar] = norm(trajectory[nn - 1][0], trajectory[nn - 1][1]);
				}
			}
			for (int i = 0; i < N; i++)
			{
				avgWalkLength += walkLength[i];
				walkLength[i] = 0; // reset for next iteration just to be sure
			}
			avgWalkLength = avgWalkLength / N;

			for (int i = 0; i < nn; i++)
			{
				trajectory[i][0] = 0;
				trajectory[i][1] = 0;
			}
			nnResults[nnIdx] = avgWalkLength;
			nnIdx++;
		}
		std::ofstream myfile;
		std::string fileName = "output/task_2/log_length_of_trj_evolution.dat";
		myfile.open(fileName);

		for (int j = 0; j < powers; j++)
			myfile << std::to_string(std::log(2) * ((long long int)j + minPowerOfTwo)) + "\t" + std::to_string(std::log(nnResults[j])) + "\n";

		myfile.close();

	/* single not returning walk
			#pragma omp parallel for
		for (int iPar = 0; iPar < N; iPar++)
		{
			for (int i = 0; i < n - 1; i++)
			{
				if (prevDir != -1)
					dirOptions[prevDir] = false;

				dir = next_direction(dirOptions);
				make_step(trajectory, dir, i, step);

				for (int q = 0; q < 4; q++)
					dirOptions[q] = true;

				prevDir = (dir + 2) % 4;
				walkLength[iPar] = norm(trajectory[n - 1][0], trajectory[n - 1][1]);
			}
		}
		for (int i = 0; i < N; i++)
			avgWalkLength += walkLength[i];
		avgWalkLength = avgWalkLength / N;

		std::cout << avgWalkLength << std::endl;
	*/
		std::cout << "Task 2 completed" << std::endl;
	}

	if (task == 3) // not crossing walk
	{
		int numOfSteps[N] = { 0 };
		double avgNumStepDeadEnd = 0;
		int numOfCollisions[N_THREAD] = { 0 };
		int collisions = 0;
		bool saveCollisionTrajectory = true; // choose whether you want to save trajectories with dead end into files

		omp_set_num_threads(N_THREAD);
		#pragma omp parallel for firstprivate(trajectory, dirOptions) 

		for (int iPar = 0; iPar < N; iPar++)
		{
			int i = 0;
			int* iPtr = &i;
			int prevDir = next_direction(dirOptions); // next_direction generated a random direction alloed by dirOptions
			int* prevDirPtr = &prevDir;	//prevDir will be updated in step_dont_cross
			bool noCollision = true;

			while (noCollision && (i < n - 1))
			{
				for (int q = 0; q < 4; q++)
					dirOptions[q] = true; // all directions are now possible
				dirOptions[prevDir] = false; // we already disallow going immediately back

				noCollision = step_dont_cross(trajectory, iPtr, step, dirOptions, prevDirPtr); 
				i++;
			}
			if (!noCollision)
			{
				numOfSteps[iPar] = i - 1;
				numOfCollisions[omp_get_thread_num()]++;

				if (saveCollisionTrajectory)
				{
					std::ofstream myfile;
					std::string fileName = "output/task_3/dead_end_" + std::to_string(iPar) + "_trj.dat";

					if (true) // (true) (i < 9) (i > 300)
					{
						myfile.open(fileName);

						for (int j = 0; j < i; j++)
						{
							myfile << std::to_string(trajectory[j][0]) + "\t" + std::to_string(trajectory[j][1]) + "\n";
						}
					}
					myfile.close();
				}
			}
		}
		for (int i = 0; i < N_THREAD; i++)
			collisions += numOfCollisions[i];

		for (int i = 0; i < N; i++)
			avgNumStepDeadEnd += numOfSteps[i];
		avgNumStepDeadEnd = avgNumStepDeadEnd / collisions;

		std::cout << "The average number of steps for non-crossing walk is: " << avgNumStepDeadEnd << std::endl;
	}
}

int next_direction(bool dirOpts[4]) // when 4 possible directions then we want to generate a random integer in the range [0, 3], if 3 possible directions, then the range is [0, 2] and for 2 possible directions we generate an integer from [0, 1]
{
	int dirLen = 0;
	for (int i = 0; i < 4; i++)
		if (dirOpts[i])
			dirLen++;

	if (dirLen == 4)
		return uni4(rng);

	else
	{
		int sortedOpts[4]; // we sort the directions such that possible ones are at the beggining of this list
		int j = 0;
		int k = 0;
		while (j < 4)
		{
			if (dirOpts[j])
				sortedOpts[j - k] = j;
			else
			{
				sortedOpts[3 - k] = j;
				k++;
			}
			j++;
		}

		if (dirLen == 3)
		{
			int idx = uni3(rng);
			return sortedOpts[idx];
		}
		else if (dirLen == 2)
		{
			int idx = uni2(rng);
			return sortedOpts[idx];
		}
		else
			return sortedOpts[0];
	}

}

void make_step(int trj[n][2], int dir, int i, int step) // makes step in the given direction
{
	switch (dir)
	{
	case 0:
		*(*(trj + i + 1) + 0) = *(*(trj + i) + 0);
		*(*(trj + i + 1) + 1) = *(*(trj + i) + 1) + step;
		break;
	case 1:
		*(*(trj + i + 1) + 0) = *(*(trj + i) + 0) + step;
		*(*(trj + i + 1) + 1) = *(*(trj + i) + 1);
		break;
	case 2:
		*(*(trj + i + 1) + 0) = *(*(trj + i) + 0);
		*(*(trj + i + 1) + 1) = *(*(trj + i) + 1) - step;
		break;
	case 3:
		*(*(trj + i + 1) + 0) = *(*(trj + i) + 0) - step;
		*(*(trj + i + 1) + 1) = *(*(trj + i) + 1);
		break;
	}
}

bool step_dont_cross(int trj[n][2], int* iPtr, int step, bool dirOpts[4], int* prevDirPtr)
{
	int dir = 0;
	bool passed = true;
	while (is_not_dead_end(dirOpts) && passed)
	{
		passed = false; // need to exit this while loop if we find the next viable step
		dir = next_direction(dirOpts);
		make_step(trj, dir, *iPtr, step);

		for (int k = 0; k < *iPtr + 1; k++)
		{
			if ((*(*(trj + k) + 0) == *(*(trj + *iPtr + 1) + 0)) && (*(*(trj + k) + 1) == *(*(trj + *iPtr + 1) + 1))) // finds out if we have already been on this position, if true than we generate a different direction, if there is no longer any possible direction then we have came to the dead end
			{
				*(dirOpts + dir) = false;
				passed = true;
				break;
			}
		}
	}
	*prevDirPtr = (dir + 2) % 4;
	return is_not_dead_end(dirOpts); // if we come to dead end, we return false and we stop the walk
}

bool is_not_dead_end(bool dirOpts[4])
{
	return (*(dirOpts + 0) || *(dirOpts + 1) || *(dirOpts + 2) || *(dirOpts + 3)); // if there is any possible direction, return true
}

double norm(int x, int y) // Euclidean norm in 2D
{
	return sqrt(x*x + y*y);
}