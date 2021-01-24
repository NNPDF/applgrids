#include <appl_grid/appl_grid.h>

#include <vector>

int main(int argc, char* argv[])
{
    // 1 = don't combine bins
    std::vector<int> bins(29, 1);
    // 2 = combine the next two bins
    bins.at(27) = 2;
    // 0 = last bin is already combined
    bins.at(28) = 0;

    appl::grid g(argv[1]);
    g.combine(bins);
    g.Write(argv[2]);
}
