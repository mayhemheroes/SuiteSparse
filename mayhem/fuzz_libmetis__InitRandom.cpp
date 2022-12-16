#include <stdint.h>
#include <stdio.h>
#include <climits>

#include <fuzzer/FuzzedDataProvider.h>

extern "C" void libmetis__InitRandom(uint64_t seed);

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    FuzzedDataProvider provider(data, size);
    int seed = provider.ConsumeIntegral<uint64_t>();
    libmetis__InitRandom(seed);

    return 0;
}