#include <string>
#include <fstream>
#include <filesystem>
#include <gtest/gtest.h>
#define private public
#include <kmtricks/kmer.hpp>
#include <kmtricks/alphabet.hpp>
#include <kmtricks/protein.hpp>
#include <kmtricks/io/kmer_file.hpp>
#include <kmtricks/io/hash_file.hpp>
#include <kmtricks/utils.hpp>

using namespace km;

// Helper function to generate random protein sequence
std::string random_protein_seq(size_t len)
{
  const char amino_acids[] = "ACDEFGHIKLMNPQRSTVWY";
  std::string result;
  result.reserve(len);
  for (size_t i = 0; i < len; ++i)
    result += amino_acids[rand() % 20];
  return result;
}

class ProteinTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    // Set protein alphabet for all tests
    Kmer<32>::set_alphabet(Alphabet::PROTEIN);
    Kmer<64>::set_alphabet(Alphabet::PROTEIN);
    Kmer<128>::set_alphabet(Alphabet::PROTEIN);
  }

  void TearDown() override
  {
    // Reset to DNA alphabet after tests
    Kmer<32>::set_alphabet(Alphabet::DNA);
    Kmer<64>::set_alphabet(Alphabet::DNA);
    Kmer<128>::set_alphabet(Alphabet::DNA);
  }
};

TEST_F(ProteinTest, alphabet_configuration)
{
  EXPECT_EQ(Kmer<32>::alphabet(), Alphabet::PROTEIN);
  EXPECT_EQ(Kmer<32>::bits_per_symbol(), 5);
  EXPECT_EQ(alphabet_bits(Alphabet::PROTEIN), 5);
  EXPECT_EQ(alphabet_bits(Alphabet::DNA), 2);
  EXPECT_EQ(alphabet_to_string(Alphabet::PROTEIN), "protein");
  EXPECT_EQ(alphabet_from_string("protein"), Alphabet::PROTEIN);
}

TEST_F(ProteinTest, kmer_encoding_decoding)
{
  // Test known protein sequences
  std::string seq1 = "ACDEFGHIKLMNPQRSTVWY";  // All 20 standard amino acids
  std::string seq2 = "AAAAACCCCCDDDDD";        // Repeated amino acids
  std::string seq3 = "MVHLTPEEKSAVTAL";        // Real protein sequence fragment

  Kmer<32> kmer1(seq1);
  EXPECT_EQ(seq1, kmer1.to_string());

  Kmer<32> kmer2(seq2);
  EXPECT_EQ(seq2, kmer2.to_string());

  Kmer<32> kmer3(seq3);
  EXPECT_EQ(seq3, kmer3.to_string());
}

TEST_F(ProteinTest, kmer_random_sequences)
{
  // Test random protein sequences of different lengths
  std::string a = random_protein_seq(12);
  std::string b = random_protein_seq(20);
  std::string c = random_protein_seq(25);

  Kmer<32> kmer_a(a);
  EXPECT_EQ(a, kmer_a.to_string());

  Kmer<64> kmer_b(b);
  EXPECT_EQ(b, kmer_b.to_string());

  Kmer<128> kmer_c(c);
  EXPECT_EQ(c, kmer_c.to_string());
}

TEST_F(ProteinTest, kmer_at_accessor)
{
  std::string seq = "ACDEFGHIKLMNPQRSTVWY";
  Kmer<32> kmer(seq);

  for (size_t i = 0; i < seq.size(); i++)
    EXPECT_EQ(seq[i], kmer.at(i));
}

TEST_F(ProteinTest, rev_comp_is_noop)
{
  // For protein kmers, reverse complement should be a no-op or identity
  std::string seq = "ACDEFGHIKLMNPQRSTVWY";
  Kmer<32> kmer(seq);
  Kmer<32> rev = kmer.rev_comp();

  // Protein reverse complement should return the same kmer
  EXPECT_EQ(kmer, rev);
}

TEST_F(ProteinTest, canonical_is_identity)
{
  // For protein kmers, canonical should return the same kmer
  std::string seq = "ACDEFGHIKLMNPQRSTVWY";
  Kmer<32> kmer(seq);
  Kmer<32> canonical = kmer.canonical();

  EXPECT_EQ(kmer, canonical);
  EXPECT_EQ(seq, canonical.to_string());
}

TEST_F(ProteinTest, kmer_comparison_operators)
{
  std::string seq_a = "AAAAAAACCCCCCC";
  std::string seq_b = "AAAAAAACCCCCCT";

  Kmer<32> kmer_a(seq_a);
  Kmer<32> kmer_b(seq_b);

  EXPECT_TRUE(kmer_a < kmer_b);
  EXPECT_FALSE(kmer_a > kmer_b);
  EXPECT_FALSE(kmer_a == kmer_b);
  EXPECT_TRUE(kmer_a != kmer_b);
  EXPECT_TRUE(kmer_a == kmer_a);
}

TEST_F(ProteinTest, kmer_slots_calculation)
{
  // Protein uses 5 bits per symbol
  // For k=12: 12 * 5 = 60 bits = 1 uint64_t
  EXPECT_EQ(kmer_slots(12, Alphabet::PROTEIN), 1);

  // For k=13: 13 * 5 = 65 bits = 2 uint64_t
  EXPECT_EQ(kmer_slots(13, Alphabet::PROTEIN), 2);

  // For k=25: 25 * 5 = 125 bits = 2 uint64_t
  EXPECT_EQ(kmer_slots(25, Alphabet::PROTEIN), 2);

  // For k=26: 26 * 5 = 130 bits = 3 uint64_t
  EXPECT_EQ(kmer_slots(26, Alphabet::PROTEIN), 3);
}

TEST_F(ProteinTest, amino_acid_encoding)
{
  // Test individual amino acid encoding
  EXPECT_EQ(protein_to_code('A'), 0);
  EXPECT_EQ(protein_to_code('C'), 1);
  EXPECT_EQ(protein_to_code('D'), 2);
  EXPECT_EQ(protein_to_code('E'), 3);
  EXPECT_EQ(protein_to_code('Y'), 19);

  // Test case insensitivity
  EXPECT_EQ(protein_to_code('a'), 0);
  EXPECT_EQ(protein_to_code('c'), 1);
  EXPECT_EQ(protein_to_code('y'), 19);

  // Test unknown/ambiguous amino acids map to X (20)
  EXPECT_EQ(protein_to_code('X'), 20);
  EXPECT_EQ(protein_to_code('B'), 20);
  EXPECT_EQ(protein_to_code('Z'), 20);
}

TEST_F(ProteinTest, io_kmer_file_round_trip)
{
  // Create temporary directory for test files
  std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "kmtricks_protein_test";
  std::filesystem::create_directories(temp_dir);
  std::string test_file = (temp_dir / "test_protein.kmer").string();

  try {
    // Write protein kmers to file
    const uint32_t k = 15;
    std::vector<std::string> sequences = {
      "ACDEFGHIKLMNPQR",
      "MVHLTPEEKSAVTAL",
      "AAAAAAACCCCCCC",
      "WWWWWWWYYYYYYY"
    };
    std::vector<uint8_t> counts = {5, 10, 15, 20};

    {
      KmerWriter<8192> writer(test_file, k, 1, 0, 0, false, Alphabet::PROTEIN);
      for (size_t i = 0; i < sequences.size(); ++i)
      {
        Kmer<32> kmer(sequences[i]);
        writer.write<32, 8>(kmer, counts[i]);
      }
    }

    // Read protein kmers from file
    {
      KmerReader<8192> reader(test_file);
      EXPECT_EQ(reader.infos().kmer_size, k);
      EXPECT_EQ(reader.infos().alphabet_id, static_cast<uint8_t>(Alphabet::PROTEIN));
      EXPECT_EQ(reader.infos().alphabet_bits_per_symbol, 5);

      size_t idx = 0;
      Kmer<32> kmer;
      kmer.set_k(k);
      uint8_t count;

      while (reader.read<32, 8>(kmer, count))
      {
        ASSERT_LT(idx, sequences.size());
        EXPECT_EQ(kmer.to_string(), sequences[idx]);
        EXPECT_EQ(count, counts[idx]);
        ++idx;
      }

      EXPECT_EQ(idx, sequences.size());
    }

    // Clean up
    std::filesystem::remove_all(temp_dir);
  }
  catch (...)
  {
    std::filesystem::remove_all(temp_dir);
    throw;
  }
}

TEST_F(ProteinTest, io_hash_file_round_trip)
{
  // Create temporary directory for test files
  std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "kmtricks_protein_test_hash";
  std::filesystem::create_directories(temp_dir);
  std::string test_file = (temp_dir / "test_protein.hash").string();

  try {
    // Write protein hashes to file
    std::vector<uint64_t> hashes = {100, 200, 300, 400, 500};
    std::vector<uint8_t> counts = {5, 10, 15, 20, 25};

    {
      HashWriter<8, 32768> writer(test_file, 1, 0, 0, false, Alphabet::PROTEIN);
      for (size_t i = 0; i < hashes.size(); ++i)
      {
        writer.write(hashes[i], counts[i]);
      }
      writer.flush();
    }

    // Read protein hashes from file
    {
      HashReader<8, 32768> reader(test_file);
      EXPECT_EQ(reader.infos().alphabet_id, static_cast<uint8_t>(Alphabet::PROTEIN));
      EXPECT_EQ(reader.infos().alphabet_bits_per_symbol, 5);

      size_t idx = 0;
      uint64_t hash;
      uint8_t count;

      while (reader.read(hash, count))
      {
        ASSERT_LT(idx, hashes.size());
        EXPECT_EQ(hash, hashes[idx]);
        EXPECT_EQ(count, counts[idx]);
        ++idx;
      }

      EXPECT_EQ(idx, hashes.size());
    }

    // Clean up
    std::filesystem::remove_all(temp_dir);
  }
  catch (...)
  {
    std::filesystem::remove_all(temp_dir);
    throw;
  }
}

TEST_F(ProteinTest, alphabet_mismatch_detection)
{
  // Create temporary directory for test files
  std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "kmtricks_protein_test_mismatch";
  std::filesystem::create_directories(temp_dir);
  std::string protein_file = (temp_dir / "protein.kmer").string();
  std::string dna_file = (temp_dir / "dna.kmer").string();

  try {
    // Write a protein kmer file
    {
      Kmer<32>::set_alphabet(Alphabet::PROTEIN);
      KmerWriter<8192> writer(protein_file, 15, 1, 0, 0, false, Alphabet::PROTEIN);
      Kmer<32> kmer("ACDEFGHIKLMNPQR");
      writer.write<32, 8>(kmer, static_cast<uint8_t>(5));
    }

    // Write a DNA kmer file
    {
      Kmer<32>::set_alphabet(Alphabet::DNA);
      KmerWriter<8192> writer(dna_file, 15, 1, 0, 0, false, Alphabet::DNA);
      Kmer<32> kmer("ACGTACGTACGTACG");
      writer.write<32, 8>(kmer, static_cast<uint8_t>(5));
    }

    // Try to aggregate files with different alphabets - should throw
    std::vector<std::string> mixed_paths = {protein_file, dna_file};
    EXPECT_THROW({
      KmerFilesAggregator<32, 8> aggregator(mixed_paths, 15);
    }, PipelineError);

    // Clean up
    std::filesystem::remove_all(temp_dir);
  }
  catch (...)
  {
    std::filesystem::remove_all(temp_dir);
    throw;
  }
}
