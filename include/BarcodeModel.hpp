#ifndef __BARCODE_MODEL_HPP__
#define __BARCODE_MODEL_HPP__

#include <iostream>
#include <string>
#include <vector>
#include <atomic>
#include <cmath>

#include "AlevinUtils.hpp"
#include "BarcodeGroup.hpp"

namespace alevin{
  namespace model{

    /*
      Caluclates Probability argmax_A {P(b | A,d)}
    */
    bool calculateAlnProbability(const std::string& s1,
                                 const std::string& s2,
                                 double& probability){
      int32_t l1{static_cast<int32_t>(s1.size())}, l2 {static_cast<int32_t>(s2.size())};

      // Shouldn't happen i.e. length diff are greater than 1
      if ( std::abs(l1-l2) > 1 ){
        probability = 0;
        return false;
      }

      char edit{'M'};
      uint32_t distance{0};
      if(std::abs(l1-l2) == 1){
        distance = 1;
      }

      for (int32_t i=0, j=0; i<l1 and j<l2; i++, j++){
        if (s1[i] != s2[j]){
          if (distance == 1){
            edit = 'F';
            break;
          }
          edit = 'S';
          distance += 1;

          if(l1-i==1 and l2-j==1){
            break;
          }
          else{
            bool flag = true;
            if(s1[i+1] == s2[j+1]){
              for(int32_t l=i+2, m=j+2; l<l1 and m<l2; l++,m++){
                if (s1[l] != s2[m]){
                  flag = false;
                  break;
                }
              }
              if (flag){
                break;
              }
            }
          }

          if(s1[i+1] == s2[j]){
            edit = 'D';
            bool flag = true;
            for(int32_t l=i+2, m=j+1; l<l1 and m<l2; l++,m++){
              if (s1[l] != s2[m]){
                flag = false;
                break;
              }
            }
            if (flag){
              break;
            }
          }

          if (s1[i] == s2[j+1]){
            edit = 'I';
            bool flag = true;
            for(int32_t l=i+1, m=j+2; l<l1 and m<l2; l++,m++){
              if (s1[l] != s2[m]){
                flag = false;
                break;
              }
            }
            if (flag){
              break;
            }
          }

          edit = 'F';
          break;
        }
      }

      switch(edit){
      case 'M': probability = 1.0;
        break;
      case 'S': probability = 0.6;//SNP
        break;
      case 'I': probability = 0.3;//Insertion
        break;
      case 'D': probability = 0.1;//deletion
        break;
      case 'F': probability = 0.0;
        return false;
      }
      return true;
    }

    template <typename T>
    std::vector<size_t> sort_indexes(const std::vector<T> &v,
                                     const std::vector<std::string>& trueBarcodes) {

      // initialize original index locations
      std::vector<size_t> idx(v.size());
      iota(idx.begin(), idx.end(), 0);

      // sort indexes based on comparing values in v
      sort(idx.begin(), idx.end(),
           [&v, &trueBarcodes](size_t i1, size_t i2) {
             if (v[i1] > v[i2])
               return true;
             else if (v[i1] < v[i2])
               return false;
             else
               if (trueBarcodes[i1] > trueBarcodes[i2])
                 return true;
               else
                 return false;
           });

      //if(v.size()>1 and v[idx[0]] == v[idx[1]]){
      //  size_t maxIdx = idx[0];
      //  size_t maxBase = 0;
      //  for (size_t qbase=1; qbase<v.size(); qbase++){
      //    size_t qIdx = idx[qbase];
      //    if (v[qIdx] < v[maxIdx]){
      //      break;
      //    }
      //    else if (v[qIdx] > v[maxIdx]) {
      //      std::cout << "ERROR in sorting the barcode probability" <<std::flush;
      //      exit(1);
      //    }
      //    if (trueBarcodes[qIdx].compare(trueBarcodes[maxIdx]) < 0){
      //      maxIdx = qIdx;
      //      maxBase = qbase;
      //    }
      //  }
      //  if (maxIdx != idx[0]){
      //    std::iter_swap(idx.begin(), idx.begin()+maxBase);
      //  }
      //}

      return idx;
    }

    void coinTossBarcodeModel(std::string barcode,
                              const std::vector<std::string>& trueBarcodes,
                              MapT& dumpPair){
      if(trueBarcodes.size() == 1){
        dumpPair.push_back(std::make_pair(trueBarcodes.front(), 1.0));
      }
      else{
        std::vector<double> probabilityVec;
        double alnProbability, probabilityNorm;

        for(const std::string trueBarcode: trueBarcodes){
          //save the sequence of the true barcodes for dumping
          bool isOneEdit = calculateAlnProbability(trueBarcode,
                                                   barcode,
                                                   alnProbability);
          if(!isOneEdit){
            std::cerr << "Barcode model receveived wrong barcode mapping."
                      << "should not happen, Please report the bug\n";
            std::cerr << trueBarcode << "\n" << barcode <<std::endl;
            std::exit(1);
          }

          probabilityVec.push_back(alnProbability);
        }//end-for

        //normalizer for the probability
        probabilityNorm = std::accumulate(probabilityVec.begin(),
                                          probabilityVec.end(),
                                          0.0);
        std::transform( probabilityVec.begin(), probabilityVec.end(),
                        probabilityVec.begin(),
                        [probabilityNorm](double n){ return n/probabilityNorm;} );

        for (auto i: sort_indexes(probabilityVec, trueBarcodes)) {
          dumpPair.push_back(std::make_pair(trueBarcodes[i], probabilityVec[i]));
        }
        //std::vector<double> cumProbVec(probabilityVec.size());
        //std::partial_sum(probabilityVec.begin(), probabilityVec.end(),
        //                 cumProbVec.begin());

        //dumping the tuple (true barcode seq, cum. prob.).
        //for(size_t i=0; i < trueBarcodes.size(); i++){
        //  dumpPair.push_back(std::make_pair(trueBarcodes[i], cumProbVec[i]));
        //}
      }//end-else
    }
  }
}

#endif
