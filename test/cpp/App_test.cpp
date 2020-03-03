#include "gtest/gtest.h"

#include <fstream>

#include "lib/namespace.hpp"

#include "app/CBase.hpp"
#include "app/CEnc.hpp"
#include "app/CBatchEnc.hpp"
#include "app/CircuitZKPVerifier.hpp"
#include "app/CircuitZKPProver.hpp"
#include "lib/math/IntegerImpl.hpp"
#include "lib/math/Matrix.hpp"
#include "lib/paillier/PaillierEncryption.hpp"
#include "lib/utils/Timer.hpp"

#include "app/App.hpp"

using namespace cryptoplus;
using namespace polyu;

namespace
{

TEST(App, Run)
{
  string filename = "../../test6.csv";
  ifstream ifile(filename);
  ofstream fs;
  if (ifile)
  {
    // ifile.cloes();
    fs.open(filename, std::ios_base::app);
  }
  else
  {
    // ifile.cloes();
    fs.open(filename);
    fs << "message count,";
    fs << "message/batch,";
    fs << "batch count,";
    fs << "byte length(bit),";
    fs << "byte length(byte),";
    fs << "message size(byte),";
    fs << "slot size(byte),";
    fs << "slot/message,";
    fs << "range proof count,";
    fs << "cipher size,";
    fs << "proof size,";
    fs << "proof size/message,";

    fs << "single encrypt circuit's N,";
    fs << "single encrypt circuit's Q,";
    fs << "batch encrypt circuit's N,";
    fs << "batch encrypt circuit's Q,";
    fs << "batch encrypt circuit's matrix m,";
    fs << "batch encrypt circuit's matrix n,";

    fs << "encryption time,";
    fs << "circuit create time,";
    fs << "value assign time,";
    fs << "commit time,";
    fs << "prove time,";
    fs << "verify time" << endl;
  }

  size_t byteLength, msgCount, rangeProofCount, slotSize, msgPerBatch;

  byteLength = 8;
  msgCount = 10;
  rangeProofCount = 10;
  slotSize = 4;
  msgPerBatch = 15;

  // vector<size_t> bls({8, 16, 32, 64, 128, 256});
  // vector<size_t> ms({10, 20, 50, 100, 200});
  vector<size_t> bls({8});
  vector<size_t> ms({10, 20, 50});
  for (size_t i = 0; i < bls.size(); i++)
  {
    byteLength = bls[i];
    auto crypto = PaillierEncryption::generate(byteLength);

    for (size_t j = 0; j < ms.size(); j++)
    {
      msgCount = ms[j];
      polyu::run(crypto, msgCount, rangeProofCount, slotSize, msgPerBatch, fs);
    }
  }

  fs.close();
}

} // namespace
