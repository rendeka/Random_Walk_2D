#define N 5000	//number of random walks
#define n 10000000	//length of each walk !WARNING! for task_1 and task_2 must be satisfied the condition: n > 2^maxPowerOftwo (Note: 2^15 = 32768)
#define N_THREAD 8		//number of threads 
#define minPowerOfTwo 2 //  minPowerOfTwo > 0 must be satisfied
#define maxPowerOfTwo 17 // in task_1 and task_2 we make a walk with length = 2^x, where x belongs to the set {minPowerOfTwo, minPowerOfTwo + 1, ..., maxPowerOfTwo - 1, maxPowerOfTwo}
#define powers maxPowerOfTwo - minPowerOfTwo

int task = 2;	//task == 1 simple walk				I have used: n = 132000 -> corresponds to maxPowerOfTwo = 17, N = 1600
				//task == 2 not returning walk		I have used: n = 132000 -> corresponds to maxPowerOfTwo = 17, N = 1600
				//task == 3 not crossing walk		I have used: n = 100000 but walks longer than 600 were never realized, N = 40000 
