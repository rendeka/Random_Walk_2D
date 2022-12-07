#pragma once
// somewhere in the file must already be #include <random>

//generating random integer
std::random_device rd; // obtain a random number from hardware
std::mt19937 rng(rd()); // seed the generator
std::uniform_int_distribution<> uni4(0, 3); // define the range
std::uniform_int_distribution<> uni3(0, 2); // define the range
std::uniform_int_distribution<> uni2(0, 1); // define the range
//end generating random integer