#include "gtest/gtest.h"

#include "Utils.hpp"
#include "EC.hpp"
#include "KeyPair.hpp"
#include "lib/utils/ECSignature.hpp"

using namespace cryptoplus;
using namespace std;

namespace
{
TEST(ECSignature, verifySignature)
{
  auto curve = EC::SECP256K1();
  string message = "foo";
  auto messageBin = Utils::stringToBinary(message);
  vector<uint8_t> seed = {0x48, 0x41, 0x48, 0x41, 0x48, 0x41};
  auto key = KeyPair::createWithSeed(curve, seed);
  vector<uint8_t> s1 = ECSignature::sign(curve, key->getPrivateKey(), messageBin);
  EXPECT_TRUE(ECSignature::verify(curve, key->getPublicKey(), messageBin, s1));
}
} // namespace
