#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <math.h>       /* pow */
#include <map>
#include <set>
#include <string>
#include <vector>
#include <queue>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "cascade.h"

class CompareDist
{
public:
    bool operator()(pair<int,double> n1,pair<int,double> n2) {
        return n1.second<n2.second;
    }
};

