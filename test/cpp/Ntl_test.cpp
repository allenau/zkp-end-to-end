#include "gtest/gtest.h"

#include <fstream>

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

#include "app/namespace.hpp"

#include "app/math/Matrix.hpp"
#include "app/utils/Timer.hpp"

#include "app/utils/ConvertUtils.hpp"
#include "app/math/MathUtils.hpp"

#include "app/App.hpp"

namespace
{

ZZ_p time_exp(const Vec<ZZ_p> &mi, const ZZ_p &r, ZZ GP_Q, ZZ_p g, Vec<ZZ_p> gi, string name, double &single_exp, uint &repeat_num) {
    
    ZZ_pPush push(GP_Q);

    ZZ_p gx;
    // Timer::start("a single exp (ns)");
    // ZZ_p ret = power(g, conv<ZZ>(r));
    // single_exp += Timer::end3("a single exp (ns)");

    Timer::start(name);
  for (size_t i = 0; i < mi.length(); i++) {
    power(gx, gi[i], conv<ZZ>(mi[i]));
    // mul(ret, ret, gx);
  }
    repeat_num = mi.length();
    single_exp = Timer::end(name, true) / repeat_num;
    //cout << "length is " << mi.length() << endl;
    
    return gx;
}

void time_mul (ZZ_p one, ZZ_p two, ZZ GP_Q, string name, double & single_mul) {
    ZZ_pPush push(GP_Q);

    ZZ_p result;
    Timer::start(name);
    for (size_t i = 0 ; i < 1000; i++) {
        mul(result, one, two);
    }
    // single_mul += Timer::end3(name);
    single_mul = Timer::end(name, true) / 1000;
    cout << "length is " << 1000 << endl;    
}

TEST(Ntl, property)
{
  cout << " === Ntl Test === " << endl;
  unsigned char bytearray[256];

  ZZ num = ConvertUtils::hexToZZ("ACB1"); //
  ZZ Y = ConvertUtils::hexToZZ("0F");
  ZZ Z;
  ZZ flag = conv<ZZ>("255");
  flag = Y;
  ZZ X = conv<ZZ>("1854952506577119873857698280485327102978808098786911742041191626187809479620209526643318380828471870727407446386445347391076366356683037156540282749078094");
  NTL::bit_and(Y, X, flag);
  NTL::bit_and(Z, X >> 4 , flag);
  NTL::bit_and(num, X >> 8 , flag);
  int x = conv<int>(Y);
  int y = conv<int>(Z);
  int z = conv<int>(num);
  cout << "0:" << x << endl;
  cout << "1:" << y << endl;
  cout << "2:" << z << endl;
  cout << "size of X:" << NumBytes(X) << endl;
  NTL::BytesFromZZ(bytearray,X, NumBytes(X));
  for (int i=0;i<256;i++) {
    cout << "X at i=" << i << "is:"<< (int)bytearray[i] << endl;
  }

  //ZZ num = conv<ZZ>("113");
  //ZZ flag = conv<ZZ>("15");
  //ZZ X = conv<ZZ>("1234");
  X = ConvertUtils::hexToZZ("FFFF");
  //ZZ Y = ConvertUtils::hexToZZ("0F");
  cout << "num:" << num << endl;
  cout << "X:" << X << endl;
  cout << "Y:" << Y << endl;
  NTL::bit_and(X, num , Y);
  cout << "Bit_AND of num and Y:" << X << endl;
  cout << "this is num:" << num << endl;
  unsigned long t = 0, k = 0;
  cout << "What is this?" << num.size() << endl;
  
  cout << "this is the last 4 bits of num:" << X << endl;
  t = (NTL::to_int(num) >> 4); 
  cout << " this is num right shifted by 4:" << t << endl;

}

TEST(Ntl, Speed_test)
{
  cout <<"Are you okay?" << endl;
  cout << "--------------- benchmark begins for testing the time cost in a signle exp. or mul. opeartion -------------" << endl;
    

    string filename = "./bench.csv";
    ifstream ifile(filename);
    ofstream fs;
    if (ifile){
    // ifile.cloes();
    fs.open(filename, std::ios_base::app);
    } else {
    // ifile.cloes();
    fs.open(filename);
    fs << "|N|" << ",";
    fs << "|Q|" << ",";
    fs << "single exp. time (ms)" << ",";
    fs << "single mul. time (ns)" << ",";
    fs << "repeat num" << "," << endl;
    }
    size_t bytelength;
    vector<size_t> bls ({128, 256}); // 1024-bit and 2048-bit
    size_t msgCount, rangeProofCount, slotSize, msgPerBatch;
    msgCount = 100;
    rangeProofCount =10; 
    slotSize = 4;
    msgPerBatch = 15;

  double single_exp, single_mul, single_mul_ave;

  for (size_t i = 0; i < bls.size(); i ++) {
    cout << "-------------data in bytelength " << bls[i] << "--------------" << endl;
    bytelength = bls[i];
    auto crypto = make_shared<PaillierEncryption>(bytelength);
 
    cout << "｜N｜： " << bytelength * 8 << endl;
    // extract system parameters from private keys
    auto GP_Q = crypto->getGroupQ(); // public parameter: group element Q
    auto GP_P = crypto->getGroupP(); // public parameter: group order p
    ZZ_p::init(GP_Q);
    auto GP_G = crypto->getGroupG();         // public parameter: group generator g
    auto pk = crypto->getPublicKey();        // public key
    auto sk1 = crypto->getPrivateElement1(); // private key component
    auto sk2 = crypto->getPrivateElement2(); // private key component
    // cout << "Q is " << GP_Q << endl;
    // cout << "n^2 is " << GP_P << endl;
    // cout << "generator is " << GP_G << endl;
    // cout << "public key n = pq is " << pk << endl;
    // cout << "p is " << sk1 << endl;
    // cout << "q is " << sk2 << endl;
    auto numBytes_Q = NumBytes(GP_Q);
    // cout << "numBytes of Group Q is " << numBytes_Q << endl;
    auto numBits_Q = numBytes_Q * 8;
    // cout << "numBits of Group Q is " << numBits_Q << endl;

    auto numBytes_N2 = NumBytes(GP_P);
    auto numBits_N2 = numBytes_N2 * 8;
    // cout << "numBytes of Group N^2 is " << numBytes_N2 << endl;
    // cout << "numBits of Group N^2 is " << numBits_N2 << endl;

    auto decryptor = make_shared<PaillierEncryption>(pk, sk1, sk2, GP_Q, GP_P, GP_G);
    auto proverCir = make_shared<CBatchEnc>(decryptor, msgCount, rangeProofCount, slotSize, msgPerBatch);
    auto giRequired = proverCir->estimateGeneratorsRequired();
    auto gi = decryptor->genGenerators(giRequired); // public paramters: generators gi for commitment scheme
    fs << "giRequired:" << giRequired <<"," << endl;

    Vec<ZZ> msg;
    for (size_t i = 0; i < msgCount; i++) {
      ZZ max = conv<ZZ>(2) << (proverCir->slotsPerMsg - 1);
      ZZ q;
      ZZ r;
      DivRem(q, r, conv<ZZ>(i), max);
      auto rStr = ConvertUtils::toBinaryString(r);
      ZZ m;
      for (size_t i = 0; i < rStr.size(); i++) {
        m = m << (proverCir->slotSize * 8);
        m += rStr[i] - '0';
      }
      msg.append(m);
    }
  
    proverCir->encrypt(msg);
    auto ljir1 = proverCir->calculateLjir();
    auto Lj = proverCir->calculateLj(ljir1);
    proverCir->wireUp(ljir1, Lj);
    proverCir->run(ljir1, Lj);
    auto prover = proverCir->generateProver(gi);

    auto Cm = proverCir->Cm;
    auto Cm_ = proverCir->Cm_;
    auto CRj = proverCir->CRj;

    ZZ_pPush push(GP_P);
    Vec<ZZ_p> randA;
    // Vec<ZZ_p> randB, randC, D;
    // ZZ_p randD;
    MathUtils::randVecZZ_p(prover->zkp->m, prover->zkp->GP_P, randA);
    // MathUtils::randVecZZ_p(prover->zkp->m, prover->zkp->GP_P, randB);
    // MathUtils::randVecZZ_p(prover->zkp->m, prover->zkp->GP_P, randC);
    // MathUtils::randVecZZ_p(prover->zkp->n, prover->zkp->GP_P, D);
    // randD = MathUtils::randZZ_p(prover->zkp->GP_P);

    //only calculate one row of exp. operation
    // Vec<ZZ_p> tmp;
    string name1;
    single_exp = 0;
    single_mul = 0;
    name1 = "EXP. time in bytelength " + to_string(bytelength);
    uint repeat_num = 0;
    ZZ_p::init(GP_Q);
    ZZ_p tmp = time_exp(prover->A[0], randA[0], GP_Q, prover->zkp->GP_G, gi, name1, single_exp, repeat_num); 
    cout << "single exp. in bytelength " << bytelength << " is " << single_exp << " ms "<< endl;

    string name2;
    name2 = "MUL. time in bytelength" + to_string(bytelength);
    time_mul(tmp, tmp, GP_Q, name2, single_mul);
    cout << "single mul. in bytelength " << bytelength << " is " << single_mul << " ns" << endl;

    fs << bytelength * 8 << ",";
    fs << numBits_Q << ",";
    fs << single_exp << ",";
    fs << single_mul << ",";
    fs << repeat_num << endl;

} //for

    // cout << "*********** test time consumed in a single exp. operation in the same group size of different power **************" << endl;


}
/*/
TEST(Ntl, Speed_test)
{

  const size_t groupTry = 5;
  const size_t runTestMax = 5000000;
  int byteLengths[] = {8, 16, 32, 64, 128, 256};
  for (size_t i = 0; i < 5; i++)
  {
    int byteLength = byteLengths[i];
    double mulThroughput = 0;
    double powThroughput = 0;
    size_t qSize = 0;
    size_t pSize = 0;
    size_t mSize = 0;

    for (size_t j = 0; j < groupTry; j++)
    {
      auto crypto = make_shared<PaillierEncryption>(byteLength);
      auto GP_Q = crypto->getGroupQ();
      auto GP_P = crypto->getGroupP();
      auto GP_G = crypto->getGroupG();
      auto pk = crypto->getPublicKey();
      auto sk1 = crypto->getPrivateElement1();
      auto sk2 = crypto->getPrivateElement2();
      qSize = max(qSize, GP_Q->toBinary().size() * 8);
      pSize = max(pSize, GP_P->toBinary().size() * 8);
      mSize = max(mSize, pk->toBinary().size() * 8);

      auto a = Random::genInteger(byteLength)->mod(GP_Q);
      auto b = Random::genInteger(byteLength)->mod(GP_Q);
      double diff, throughput;

      ZZ zzQ = conv<ZZ>(GP_Q->toString().c_str());
      ZZ_p::init(zzQ);
      ZZ_p zzG = conv<ZZ_p>(GP_G->toString().c_str());
      ZZ_p zzA = conv<ZZ_p>(a->toString().c_str());
      ZZ_p zzB = conv<ZZ_p>(b->toString().c_str());

      size_t runTest = runTestMax;
      Timer::start("mod");
      for (size_t i = 0; i < runTestMax; i++)
      {
        zzA = zzA * zzB;
      }
      diff = Timer::end("mod", true);
      throughput = runTest / diff;
      mulThroughput += throughput;
      EXPECT_TRUE(a->gt(Integer::ZERO()));

      zzB = zzA;
      zzA = zzG;
      runTest = runTestMax / byteLength / 100;
      Timer::start("pow");
      ZZ zzb = conv<ZZ>(zzB);
      for (size_t i = 0; i < runTest; i++)
      {
        zzA = power(zzA, zzb);
      }
      diff = Timer::end("pow", true);
      throughput = runTest / diff;
      powThroughput += throughput;
      EXPECT_TRUE(a->gt(Integer::ZERO()));
    }

    cout << "key length: " << mSize << " bit" << endl;
    cout << "group p size: " << pSize << " bit" << endl;
    cout << "group Q size: " << qSize << " bit" << endl;
    cout << "average throughput: " << (mulThroughput / groupTry) << " mul/sec" << endl;
    cout << "average throughput: " << (powThroughput / groupTry) << " pow/sec" << endl;
  }
}
//*/

} // namespace
