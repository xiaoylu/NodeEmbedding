#include "hierar_infer.h"

HierarInfer::HierarInfer(string txtfile)
{
    loadData(txtfile);
}

void HierarInfer::loadData(string txtfile)
{
    cout << "read txt file "<< txtfile <<"...." << endl;
    ifstream file(txtfile.c_str());
    string str; 
    int node;

    // Read the node names.
    // format in the txt file:
    //   node1,1
    //   node2,2
    //     ...
    //   nodeN,N
    while (getline(file, str))
    {
        if (str.empty()) break;
        vector<string> elems;
        split(str, ',', elems);
        node = atoi((elems[0]).c_str());
		node2str.insert(pair<int, string>(node, elems[1])); 
        node2casc[node] = new set<string>();
    }
    num_nodes = node2str.size(); 

    while (getline(file, str))
    {
#ifdef DEBUG_MODE
        cout << str << endl;
#endif
        if (str.empty()) continue; //skip the empty lines
        vector<string>::iterator it;
        vector<string> elems;
        split(str, ';', elems);
        vector<string> elems2;
        split(elems[1], ',', elems2);

        string casc = elems[0];
        casc_list[casc] = new Cascade(elems[0]);

        for (it = elems2.begin(); it != elems2.end(); ++it)
        {
            int node = atoi((*it).c_str());
            ++it; // go to the time

            double time = atof((*it).c_str());
            node2casc[node]->insert(casc);
            casc_list[casc]->addInfection(node, time);
        }
	}
  
    num_cascades = casc_list.size();
}

HierarInfer::~HierarInfer()
{
   // release memory 
    for(auto& it : casc_list)
    {
      Cascade * casc = it.second;
      delete casc;
    }
}

void HierarInfer::filter(double min_duration, double max_duration)
{
    cout << "A total of " << casc_list.size() << " raw cascades" << endl;
    cout << "Filtering" << endl;
    double min_t, max_t, duration;

    for (auto it = casc_list.begin(); it != casc_list.end(); /*no increment*/ )
    {
        Cascade * casc = it->second;
        duration = casc->getDuration(min_t, max_t);
        if (casc->getLength() < 2
            || duration > max_duration || duration < min_duration)
        {
           casc_list.erase(it++); // Delete illegal cascade
        }
        else ++it;
    }
    cout << "A total of " << casc_list.size() << " cascades for analysis" << endl;
}

double HierarInfer::averageDuration()
{
    double min_t, max_t, sum_duration = 0;
    for (auto& it : casc_list)
    {
        Cascade * casc = it.second;
        sum_duration += casc->getDuration(min_t, max_t);
    }
    return sum_duration / casc_list.size();
}

void HierarInfer::normalize()
{
    double ave_t = averageDuration();
    double min_t, max_t;
    for (auto& it : casc_list)
    {
        Cascade * casc = it.second;
        casc->getDuration(min_t, max_t);
#ifdef DEBUG_MODE
        cout << min_t << ";" << max_t << endl;
#endif
        casc->alignTime(min_t);
        casc->normalizeTime(ave_t);
    }

    OW = averageDuration();
    cout << "Average cascade duration " << ave_t; 
    cout << "-->" << OW << endl;

    constructCascIDs();
}

void HierarInfer::constructCascIDs()
{
    // construct the ID lists
    casc_IDs.clear();
    for (auto& pair : casc_list)
    {
        casc_IDs.push_back(pair.first);
    }
}

double HierarInfer::loglikelihood()
{
    double ret = 0;
    int tid;

    Matrix2::instance().sumB();

    #pragma omp parallel private(tid)
    {
        tid = omp_get_thread_num();
        //printf("Thread %d starting...\n",tid);
    
        int size = casc_IDs.size(); 

        #pragma omp parallel for reduction(+:ret)
        for (int i = 0;i < size;++i)
        {
            Cascade * casc = casc_list[casc_IDs[i]];
            //cout << "Thread " << tid << ", log-likelihood=";
            //cout  << casc->logLikelihood() << endl;
            ret = ret + casc->logLikelihood();
        }
    }

    return ret;
}

void HierarInfer::infer(int dimension, double OW,
                        double stepsize, double epsilon,
                        int negative_sampling, int max_iterations,
                        int early_stop)
{
/* dimension:   length of the vectors
   OW:          the observation window of a cascades, i.e. how long the system monitors a cascade. OW controls the penalization on the "negative" links along which the infection does NOT happen within the duration OW. More specifically, this design considers the ZEROs in the adjacent matrix of the propagation graph because negative example are also important.
   stepsize:    the stepsize of the gradient descent
   epsilon:     the minimum value 
   negative_sampling: the number of negative samples chosen
   max_iterations : the maximum number of iterations
*/

    //==============sampling cascades===============
    cout << "Memory Allocation For Fast Inferrence......" << endl;
    Matrix2::instance().init(num_nodes, dimension);
    for (auto& it : casc_list)
    {
        Cascade * casc = it.second;
        casc->setParams(dimension, OW, stepsize, epsilon, negative_sampling);
        casc->createBuffer();
    }

    // gradient descent
    cout << "Learning..." << endl;
    double max_logl = loglikelihood();
    int early_stop_counter = early_stop;
    for (int round = 0;round < max_iterations; ++round) 
    {
      double logl = loglikelihood();

      if (max_logl <= logl)
      {  
        max_logl = logl;
        early_stop_counter = early_stop;
      }
      else
      {
        if (--early_stop_counter == 0)
        {
           cout << "Early stop." << endl;
           break;
        }
      }

      cout << "---------logLikelihood=" << logl << "---------" << endl;

      cout << "Round " << round << "Start" << endl;
      int size = casc_IDs.size(); 

      #pragma omp parallel for
      for (int i = 0;i < size;++i)
      {
          Cascade * casc = casc_list[casc_IDs[i]];
          casc->update();
      }
    }
    cout << "========Final LogLikelihood=" << loglikelihood() << "========" << endl;
}

void HierarInfer::output(string output_path)
{
  cout << "output the file to " << output_path << endl;
  Matrix2::instance().printMaxComponent();
  Matrix2::instance().writeMat( output_path ); 
}

//=================================
// helper function
void HierarInfer::split(const string &s, char delim, vector<string> &elems)
{
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) 
    {
        elems.push_back(item);
    }
}
