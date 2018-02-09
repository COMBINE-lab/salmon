# How to run unit test
```
./run.sh <PATH to a new Salmon Binary containg Alevin>
```
# Contents of unit_test_data.tar.gz

* *alevin*
   * __alevin.log__: logs of a run of salmon with the test data  
   * __counts.mat__: Cell(row)xGene(Column) counts  
   * __eq_classes.txt__: Global eqClass  
   * *cell*  
        * __cell__eq_classes.txt__: EqClass for one Cell  
        * __quant.sf__: Abundance of one cell  

* *src-py*
   * __get_correlation.py__: python script to get correlation of one cell
