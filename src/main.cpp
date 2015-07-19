/*
 CtuCopy 3 - Universal feature extractor and speech enhancer.
 Copyright 2012 Petr Fousek, FEE CTU Prague

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

#include "io/stdafx.h"
#include "io/opts.h"
#include "io/in.h"
#include "io/batch.h"

using namespace std;
int main (int argc, char *argv[])
{
  char version[] = VERSION;
  opts *o = NULL;
  BATCH *batch = NULL;

  // ------- INIT -------
  try{
      o = new opts(argc,argv,version);
  }
  catch (const char * p){
     cerr <<p<<endl; return -1;
  }catch (std::string p){
     cerr <<p<<endl; return -1;
  }

  try{
     batch = new BATCH(o);
  }
  catch (const char * p){
     cerr <<p<<endl;
     return -1;
  }catch (std::string p){
     cerr <<p<<endl;
     return -1;
  }

  // ------ PROCESS --------
  try {
     batch->process();
  }catch (const char * p){
     cerr << p << endl;
     return -1;
  }catch (string p){
     cerr <<p<<endl;
     return -1;
  }

  // -------- END -------
  delete batch;
  delete o;

  return 0;
}



