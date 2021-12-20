#include <cstdlib>
#include <fstream>
#include <string>

int main(int __argc, char **__argv) {
  srand(time(NULL));
  int __votex = atoi(__argv[1]);
  int __edges = atoi(__argv[2]);
  std::ofstream fos(__argv[1]);
  fos << __votex << "\n";
  fos << __edges << "\n";
  for (int i = 0; i < __edges; ++i) {
    int __from = rand() % __votex;
    int __to;
    do {
      __to = rand() % __votex;
    } while (__from == __to);
    double __distance = ((double) (rand() % 10000)) / 10000;
    fos << __from << " " << __to << " " << __distance << "\n";
  }
  fos.flush();
  fos.close();
}