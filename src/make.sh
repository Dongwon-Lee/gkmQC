g++ -Wall -Wconversion -O3 -fPIC -c *.cpp
g++ -Wall -Wconversion -O3 -fPIC *.o -lm -lpthread -o GkmKernel.so -shared