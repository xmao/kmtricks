/*****************************************************************************
 *   kmtricks
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <vector>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <gatb/gatb_core.hpp>

#include <kmtricks/alphabet.hpp>
#include <kmtricks/cmd/all.hpp>
#include <kmtricks/cmd/cmd_common.hpp>
#include <kmtricks/cmd/count.hpp>
#include <kmtricks/cmd/repart.hpp>
#include <kmtricks/cmd/superk.hpp>
#include <kmtricks/exceptions.hpp>
#include <kmtricks/hash.hpp>
#include <kmtricks/histogram.hpp>
#include <kmtricks/io/hist_file.hpp>
#include <kmtricks/io/hash_file.hpp>
#include <kmtricks/io/kff_file.hpp>
#include <kmtricks/io/kmer_file.hpp>
#include <kmtricks/io/superk_storage.hpp>
#include <kmtricks/io/vector_file.hpp>
#include <kmtricks/kmer.hpp>
#include <kmtricks/kmer_hash.hpp>
#include <kmtricks/kmdir.hpp>
#include <kmtricks/itask.hpp>
#include <kmtricks/task.hpp>
#include <kmtricks/task_pool.hpp>
#include <kmtricks/utils.hpp>

namespace km {

namespace fs = std::filesystem;

enum class ProteinOutput
{
  KMER,
  KFF,
  HASH,
  VECTOR
};

// Histogram configuration
constexpr uint32_t HIST_MIN_COUNT = 1;
constexpr uint32_t HIST_MAX_COUNT = 255;

struct ProteinMask
{
  size_t full_words {0};
  size_t clear_start {0};
  uint64_t rem_mask {0};
  bool has_rem {false};
};

inline ProteinMask make_kmer_mask(size_t kmer_size, uint8_t bits_per_symbol)
{
  ProteinMask mask;
  const size_t used_bits = kmer_size * bits_per_symbol;
  mask.full_words = used_bits / 64;
  const size_t rem_bits = used_bits % 64;
  mask.has_rem = rem_bits != 0;
  if (mask.has_rem)
    mask.rem_mask = (rem_bits == 64) ? ~0ULL : ((1ULL << rem_bits) - 1);
  mask.clear_start = mask.full_words + (mask.has_rem ? 1 : 0);
  return mask;
}

template<size_t MAX_K>
inline void apply_kmer_mask(Kmer<MAX_K>& kmer, const ProteinMask& mask)
{
  auto* data = kmer.get_data64_unsafe();
  if (mask.has_rem)
    data[mask.full_words] &= mask.rem_mask;
  for (size_t i = mask.clear_start; i < Kmer<MAX_K>::m_max_data; ++i)
    data[i] = 0;
}

template<typename Func>
inline void for_each_sequence(const std::string& files, Func&& func)
{
  IBank* bank = Bank::open(files); LOCAL(bank);
  Iterator<Sequence>* it = bank->iterator(); LOCAL(it);
  for (it->first(); !it->isDone(); it->next())
  {
    const Sequence& seq = it->item();
    const char* data = seq.getDataBuffer();
    const size_t len = seq.getDataSize();
    if (len > 0)
      func(data, len);
  }
}

inline uint32_t normalize_nb_parts(uint32_t nb_parts)
{
  if (nb_parts == 0)
    nb_parts = 4;
  if (nb_parts < 4)
    nb_parts = 4;
  return nb_parts;
}

template<size_t MAX_K>
inline Configuration write_protein_config(uint32_t kmer_size,
                                          uint32_t minim_size,
                                          uint32_t nb_parts,
                                          uint32_t nb_threads,
                                          uint32_t max_memory,
                                          uint64_t bloom_size)
{
  Configuration config;
  config._kmerSize = kmer_size;
  config._minim_size = minim_size;
  config._repartitionType = 0;
  config._minimizerType = 0;
  config._max_memory = max_memory;
  config._nbCores = nb_threads;
  config._nb_partitions = nb_parts;
  config._nb_bits_per_kmer = static_cast<uint16_t>(kmer_size * Kmer<MAX_K>::bits_per_symbol());
  config._nb_cached_items_per_core_per_part = 0;
  config._abundanceUserNb = 1;

  Storage* config_storage =
    StorageFactory(STORAGE_FILE).create(KmDir::get().m_config_storage, true, false);
  LOCAL(config_storage);
  config.save(config_storage->getGroup("gatb"));

  HashWindow hw(bloom_size, nb_parts, minim_size);
  hw.serialize(KmDir::get().m_hash_win);

  KmDir::get().init_part(nb_parts);
  return config;
}

template<size_t MAX_K, size_t MAX_C>
class ProteinCountTask : public ITask
{
  using count_type = typename selectC<MAX_C>::type;

public:
  ProteinCountTask(const std::string& sample_id,
                   uint32_t sample_idx,
                   const std::string& files,
                   uint32_t kmer_size,
                   uint32_t nb_parts,
                   uint32_t abundance_min,
                   bool lz4,
                   ProteinOutput output,
                   const std::vector<uint32_t>& partitions,
                   uint64_t hash_window,
                   hist_t hist,
                   bool clear = true)
    : ITask(3),
      m_sample_id(sample_id),
      m_sample_idx(sample_idx),
      m_files(files),
      m_kmer_size(kmer_size),
      m_nb_parts(nb_parts),
      m_ab_min(abundance_min),
      m_lz4(lz4),
      m_output(output),
      m_partitions(partitions),
      m_hash_window(hash_window),
      m_hist(hist),
      m_clear(clear)
  {
    m_part_mask.assign(m_nb_parts, false);
    for (auto part : m_partitions)
      m_part_mask[part] = true;
  }

  void preprocess() override {}
  void postprocess() override
  {
    this->m_finish = true;
    this->exec_callback();
  }

  void exec() override
  {
    spdlog::debug("[exec] - ProteinCountTask - S={}", m_sample_id);
    const bool hash_output = (m_output == ProteinOutput::HASH || m_output == ProteinOutput::VECTOR);
    if (hash_output && m_hash_window == 0)
      throw ConfigError("Hash output requested but bloom window is zero.");

    const size_t slots = Kmer<MAX_K>::slots_for_kmer(m_kmer_size);
    const size_t record_size = slots * sizeof(uint64_t);

    std::vector<std::string> temp_paths(m_nb_parts);
    std::vector<std::ofstream> temp_out(m_nb_parts);
    for (auto part : m_partitions)
    {
      temp_paths[part] = fmt::format("{}/partition_{}/{}.protein.tmp",
                                     KmDir::get().m_counts_storage,
                                     part,
                                     m_sample_id);
      temp_out[part].open(temp_paths[part], std::ios::binary | std::ios::out);
      check_fstream_good(temp_paths[part], temp_out[part]);
    }

    const uint8_t bits = Kmer<MAX_K>::bits_per_symbol();
    ProteinMask mask = make_kmer_mask(m_kmer_size, bits);
    Kmer<MAX_K> kmer;
    kmer.set_k(m_kmer_size);
    typename KmerHashers<0>::template Hasher<MAX_K> base_hasher;

    for_each_sequence(m_files, [&](const char* data, size_t len) {
      kmer.zero();
      size_t filled = 0;
      for (size_t i = 0; i < len; ++i)
      {
        const uint8_t code = char_to_code(Alphabet::PROTEIN, static_cast<unsigned char>(data[i]));
        kmer = (kmer << bits) + static_cast<uint64_t>(code);
        if (filled < m_kmer_size)
        {
          filled++;
          if (filled < m_kmer_size)
            continue;
        }
        apply_kmer_mask(kmer, mask);
        const uint64_t hash_base = base_hasher(kmer);
        const uint32_t part = static_cast<uint32_t>(hash_base % m_nb_parts);
        if (!m_part_mask[part])
          continue;
        if (hash_output)
        {
          const uint64_t hash = (hash_base % m_hash_window) + (m_hash_window * part);
          temp_out[part].write(reinterpret_cast<const char*>(&hash), sizeof(hash));
        }
        else
        {
          temp_out[part].write(reinterpret_cast<const char*>(kmer.get_data64()), record_size);
        }
      }
    });

    for (auto part : m_partitions)
      temp_out[part].close();

    for (auto part : m_partitions)
    {
      if (m_output == ProteinOutput::KMER || m_output == ProteinOutput::KFF)
        write_kmer_partition(part, temp_paths[part], record_size);
      else
        write_hash_partition(part, temp_paths[part]);
      if (m_clear)
        fs::remove(temp_paths[part]);
    }

    spdlog::debug("[done] - ProteinCountTask - S={}", m_sample_id);
  }

private:
  void write_kmer_partition(uint32_t part, const std::string& temp_path, size_t record_size)
  {
    std::vector<Kmer<MAX_K>> kmers;
    if (fs::exists(temp_path))
    {
      const size_t bytes = fs::file_size(temp_path);
      if (record_size != 0)
        kmers.reserve(bytes / record_size);
      std::ifstream in(temp_path, std::ios::binary | std::ios::in);
      check_fstream_good(temp_path, in);
      std::vector<uint64_t> buffer(record_size / sizeof(uint64_t), 0);
      while (in.read(reinterpret_cast<char*>(buffer.data()), record_size))
      {
        kmers.emplace_back();
        std::memcpy(kmers.back().get_data64_unsafe(), buffer.data(), record_size);
      }
    }

    const std::string out_path = KmDir::get().get_count_part_path(
      m_sample_id, part, m_lz4, m_output == ProteinOutput::KFF ? KM_FILE::KFF : KM_FILE::KMER);

    const count_type max_c = std::numeric_limits<count_type>::max();

    if (m_output == ProteinOutput::KFF)
    {
      kff_w_t<MAX_C> writer = std::make_shared<KffWriter<MAX_C>>(out_path,
                                                                 m_kmer_size,
                                                                 Alphabet::PROTEIN);
      if (!kmers.empty())
      {
        std::sort(kmers.begin(), kmers.end());
        size_t idx = 0;
        while (idx < kmers.size())
        {
          size_t next = idx + 1;
          uint32_t count = 1;
          while (next < kmers.size() && kmers[next] == kmers[idx])
          {
            count++;
            next++;
          }
          if (m_hist) m_hist->inc(count);
          if (count >= m_ab_min)
          {
            const count_type out_c = count > max_c ? max_c : static_cast<count_type>(count);
            writer->template write<MAX_K>(kmers[idx], out_c);
          }
          idx = next;
        }
      }
      writer->close();
      return;
    }

    kw_t<8192> writer = std::make_shared<KmerWriter<8192>>(out_path,
                                                           m_kmer_size,
                                                           requiredC<MAX_C>::value/8,
                                                           m_sample_idx,
                                                           part,
                                                           m_lz4,
                                                           Alphabet::PROTEIN);
    if (kmers.empty())
      return;

    std::sort(kmers.begin(), kmers.end());
    size_t idx = 0;
    while (idx < kmers.size())
    {
      size_t next = idx + 1;
      uint32_t count = 1;
      while (next < kmers.size() && kmers[next] == kmers[idx])
      {
        count++;
        next++;
      }
      if (m_hist) m_hist->inc(count);
      if (count >= m_ab_min)
      {
        const count_type out_c = count > max_c ? max_c : static_cast<count_type>(count);
        writer->template write<MAX_K, MAX_C>(kmers[idx], out_c);
      }
      idx = next;
    }
  }

  void write_hash_partition(uint32_t part, const std::string& temp_path)
  {
    std::vector<uint64_t> hashes;
    if (fs::exists(temp_path))
    {
      const size_t bytes = fs::file_size(temp_path);
      hashes.reserve(bytes / sizeof(uint64_t));
      std::ifstream in(temp_path, std::ios::binary | std::ios::in);
      check_fstream_good(temp_path, in);
      uint64_t hash = 0;
      while (in.read(reinterpret_cast<char*>(&hash), sizeof(hash)))
        hashes.push_back(hash);
    }

    const count_type max_c = std::numeric_limits<count_type>::max();

    if (m_output == ProteinOutput::HASH)
    {
      const std::string out_path = KmDir::get().get_count_part_path(
        m_sample_id, part, m_lz4, KM_FILE::HASH);
      hw_t<MAX_C, 32768> writer = std::make_shared<HashWriter<MAX_C, 32768>>(
        out_path, requiredC<MAX_C>::value/8, m_sample_idx, part, m_lz4, Alphabet::PROTEIN);
      if (!hashes.empty())
      {
        std::sort(hashes.begin(), hashes.end());
        size_t idx = 0;
        while (idx < hashes.size())
        {
          size_t next = idx + 1;
          uint32_t count = 1;
          while (next < hashes.size() && hashes[next] == hashes[idx])
          {
            count++;
            next++;
          }
          if (m_hist) m_hist->inc(count);
          if (count >= m_ab_min)
          {
            const count_type out_c = count > max_c ? max_c : static_cast<count_type>(count);
            writer->write(hashes[idx], out_c);
          }
          idx = next;
        }
      }
      writer->flush();
      return;
    }

    const std::string out_path = KmDir::get().get_count_part_path(
      m_sample_id, part, m_lz4, KM_FILE::VECTOR);
    bvw_t<8192> writer = std::make_shared<BitVectorWriter<8192>>(
      out_path, m_hash_window, m_sample_idx, part, m_lz4, Alphabet::PROTEIN);
    std::vector<uint8_t> vec(NBYTES(m_hash_window), 0);
    if (!hashes.empty())
    {
      std::sort(hashes.begin(), hashes.end());
      size_t idx = 0;
      while (idx < hashes.size())
      {
        size_t next = idx + 1;
        uint32_t count = 1;
        while (next < hashes.size() && hashes[next] == hashes[idx])
        {
          count++;
          next++;
        }
        if (m_hist) m_hist->inc(count);
        if (count >= m_ab_min)
          BITSET(vec, hashes[idx] - (m_hash_window * part));
        idx = next;
      }
    }
    writer->write(vec);
  }

private:
  std::string m_sample_id;
  uint32_t m_sample_idx;
  std::string m_files;
  uint32_t m_kmer_size;
  uint32_t m_nb_parts;
  uint32_t m_ab_min;
  bool m_lz4;
  bool m_clear;
  ProteinOutput m_output;
  std::vector<uint32_t> m_partitions;
  std::vector<bool> m_part_mask;
  uint64_t m_hash_window;
  hist_t m_hist;
};

template<size_t MAX_K, size_t MAX_C>
inline void run_protein_all(all_options_t opt)
{
  KmDir::get().init(opt->dir, opt->fof, true);
  opt->dump(KmDir::get().m_options);

  const uint32_t nb_parts = normalize_nb_parts(opt->nb_parts);
  Configuration config = write_protein_config<MAX_K>(opt->kmer_size,
                                                     opt->minim_size,
                                                     nb_parts,
                                                     opt->nb_threads,
                                                     opt->max_memory,
                                                     opt->bloom_size);

  opt->m_ab_min_vec.resize(KmDir::get().m_fof.size(), opt->m_ab_min);
  if (opt->restrict_to_list.empty())
  {
    if (opt->restrict_to != 1.0)
    {
      std::vector<uint32_t> parts;
      for (uint32_t i = 0; i < config._nb_partitions; ++i)
        parts.push_back(i);
      std::random_device rd;
      std::mt19937 g(rd());
      std::shuffle(parts.begin(), parts.end(), g);
      size_t n = config._nb_partitions * opt->restrict_to;
      if (n == 0) n = 1;
      for (size_t i = 0; i < n; ++i)
        opt->restrict_to_list.push_back(parts[i]);
    }
    else
    {
      for (uint32_t i = 0; i < config._nb_partitions; ++i)
        opt->restrict_to_list.push_back(i);
    }
  }
  else
  {
    for (auto part : opt->restrict_to_list)
    {
      if (part >= config._nb_partitions)
        throw PipelineError(fmt::format("Ask to process part {} but nb_partitions is {}",
                                        part, config._nb_partitions));
    }
  }

  std::vector<hist_t> hists;
  hists.resize(KmDir::get().m_fof.size(), nullptr);
  for (size_t i = 0; i < KmDir::get().m_fof.size(); ++i)
    hists[i] = opt->hist ? std::make_shared<KHist>(i, opt->kmer_size, HIST_MIN_COUNT, HIST_MAX_COUNT) : nullptr;

  const ProteinOutput output = opt->count_format == COUNT_FORMAT::HASH
    ? ProteinOutput::HASH
    : (opt->kff ? ProteinOutput::KFF : ProteinOutput::KMER);

  HashWindow hw(KmDir::get().m_hash_win);
  const uint64_t hash_window = hw.get_window_size_bits();

  if (opt->until == COMMAND::REPART || opt->until == COMMAND::SUPERK)
  {
    Eraser::get().join();
    return;
  }

  TaskPool pool(opt->nb_threads);
  for (auto id : KmDir::get().m_fof)
  {
    const std::string& sample_id = std::get<0>(id);
    const uint32_t sample_idx = KmDir::get().m_fof.get_i(sample_id);
    const std::string files = KmDir::get().m_fof.get_files(sample_id);
    const uint32_t a_min = std::get<2>(id) == 0 ? opt->c_ab_min : std::get<2>(id);
    task_t task = std::make_shared<ProteinCountTask<MAX_K, MAX_C>>(
      sample_id,
      sample_idx,
      files,
      opt->kmer_size,
      config._nb_partitions,
      a_min,
      opt->lz4,
      output,
      opt->restrict_to_list,
      hash_window,
      hists[sample_idx]);
    pool.add_task(task);
  }
  pool.join_all();

  if (opt->hist)
  {
    for (auto& h : hists)
      if (h) h->merge_clones();
    for (auto& h : hists)
      if (h) HistWriter(KmDir::get().get_hist_path(KmDir::get().m_fof.get_id(h->idx())),
                        *h,
                        false,
                        Kmer<MAX_K>::alphabet());
  }

  if (opt->until == COMMAND::COUNT || opt->kff)
  {
    Eraser::get().join();
    return;
  }

  if (opt->m_ab_float)
    opt->m_ab_min_vec = compute_merge_thresholds(hists, opt->m_ab_min_f,
                                                 KmDir::get().get_merge_th_path());

  TaskPool merge_pool(opt->nb_threads);
  for (auto part : opt->restrict_to_list)
  {
    task_t task = nullptr;
    if (opt->count_format == COUNT_FORMAT::KMER)
    {
      task = std::make_shared<KmerMergeTask<MAX_K, MAX_C>>(
        part, opt->m_ab_min_vec, config._kmerSize, opt->r_min, opt->save_if,
        opt->lz4, opt->mode, opt->format, !opt->keep_tmp);
    }
    else
    {
      task = std::make_shared<HashMergeTask<MAX_C>>(
        part, opt->m_ab_min_vec, opt->r_min, opt->save_if, opt->lz4, opt->mode,
        opt->format, hw, !opt->keep_tmp, opt->bwidth);
    }
    merge_pool.add_task(task);
  }
  merge_pool.join_all();
  Eraser::get().join();
}

template<size_t MAX_K, size_t MAX_C>
inline void run_protein_count(count_options_t opt)
{
  KmDir::get().init(opt->dir, "", false);
  Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
  LOCAL(config_storage);
  Configuration config = Configuration();
  config.load(config_storage->getGroup("gatb"));
  KmDir::get().init_part(config._nb_partitions);

  std::vector<uint32_t> partitions;
  if (opt->partition_id != -1)
  {
    if (opt->partition_id >= static_cast<int32_t>(config._nb_partitions))
      throw ConfigError(fmt::format("Ask to process partition {} but nb_partitions is {}",
                                    opt->partition_id, config._nb_partitions));
    partitions.push_back(static_cast<uint32_t>(opt->partition_id));
  }
  else
    for (uint32_t p = 0; p < config._nb_partitions; ++p) partitions.push_back(p);

  HashWindow hw(KmDir::get().m_hash_win);
  const uint64_t hash_window = hw.get_window_size_bits();

  ProteinOutput output = ProteinOutput::KMER;
  if (opt->format == "kff")
    output = ProteinOutput::KFF;
  else if (opt->format == "hash")
    output = ProteinOutput::HASH;
  else if (opt->format == "vector")
    output = ProteinOutput::VECTOR;

  const uint32_t sample_idx = KmDir::get().m_fof.get_i(opt->id);
  const std::string files = KmDir::get().m_fof.get_files(opt->id);
  hist_t hist = opt->hist ? std::make_shared<KHist>(sample_idx, config._kmerSize, HIST_MIN_COUNT, HIST_MAX_COUNT) : nullptr;

  ProteinCountTask<MAX_K, MAX_C> task(opt->id,
                                      sample_idx,
                                      files,
                                      config._kmerSize,
                                      config._nb_partitions,
                                      opt->c_ab_min,
                                      opt->lz4,
                                      output,
                                      partitions,
                                      hash_window,
                                      hist,
                                      opt->clear);
  task.preprocess();
  task.exec();
  task.postprocess();

  if (opt->hist && hist)
  {
    hist->merge_clones();
    HistWriter(KmDir::get().get_hist_path(opt->id), *hist, false, Kmer<MAX_K>::alphabet());
  }
}

template<size_t MAX_K>
inline void run_protein_repart(repart_options_t opt)
{
  KmDir::get().init(opt->dir, opt->fof, true);
  const uint32_t nb_parts = normalize_nb_parts(opt->nb_parts);
  write_protein_config<MAX_K>(opt->kmer_size,
                              opt->minim_size,
                              nb_parts,
                              opt->nb_threads,
                              0,
                              opt->bloom_size);
  spdlog::info("Protein repartition skipped (hash-based partitions).");
}

template<size_t MAX_K>
inline void run_protein_superk(superk_options_t opt)
{
  KmDir::get().init(opt->dir, "", false);

  Storage* config_storage = StorageFactory(STORAGE_FILE).load(KmDir::get().m_config_storage);
  LOCAL(config_storage);
  Configuration config = Configuration();
  config.load(config_storage->getGroup("gatb"));

  // Protein super-k-mers use hash-based minimizers (no canonical form needed)
  // Extract k-mers grouped by minimizer for memory efficiency
  Kmer<MAX_K>::set_alphabet(Alphabet::PROTEIN);

  const std::string& sample_id = opt->id;
  const std::string files = KmDir::get().m_fof.get_files(sample_id);
  const uint32_t kmer_size = config._kmerSize;
  const uint32_t minim_size = config._minim_size;

  spdlog::info("Generating protein super-k-mers for {} (k={}, m={})", sample_id, kmer_size, minim_size);

  // Open output storage for super-k-mers (partitioned by minimizer hash)
  std::unordered_set<int> all_partitions;
  for (uint32_t p = 0; p < config._nb_partitions; ++p)
    all_partitions.insert(p);

  SuperKStorageWriter* writer = new SuperKStorageWriter(
    KmDir::get().get_superk_path(sample_id), "pskp", config._nb_partitions, opt->lz4, all_partitions);

  // Process sequences and extract super-k-mers
  for_each_sequence(files, [&](const char* data, size_t len) {
    if (len < kmer_size) return;

    // Build k-mers and track minimizers
    std::vector<Kmer<MAX_K>> kmers;
    std::vector<uint64_t> minimizers;

    Kmer<MAX_K> kmer;
    kmer.set_k(kmer_size);
    const uint8_t bits = Kmer<MAX_K>::bits_per_symbol();
    ProteinMask mask = make_kmer_mask(kmer_size, bits);

    kmer.zero();
    size_t filled = 0;

    for (size_t i = 0; i < len; ++i)
    {
      const uint8_t code = char_to_code(Alphabet::PROTEIN, static_cast<unsigned char>(data[i]));
      kmer = (kmer << bits) + static_cast<uint64_t>(code);

      if (filled < kmer_size)
      {
        filled++;
        if (filled < kmer_size) continue;
      }

      apply_kmer_mask(kmer, mask);
      Mmer minim = kmer.minimizer(minim_size);

      kmers.push_back(kmer);
      minimizers.push_back(minim.value());
    }

    // Group k-mers by consecutive minimizers to form super-k-mers
    if (kmers.empty()) return;

    size_t start = 0;
    for (size_t i = 1; i <= kmers.size(); ++i)
    {
      // Check if minimizer changed or reached end
      if (i == kmers.size() || minimizers[i] != minimizers[start])
      {
        // Write super-k-mer for range [start, i)
        const uint32_t part = static_cast<uint32_t>(minimizers[start] % config._nb_partitions);
        const size_t superk_len = i - start + kmer_size - 1;
        const uint8_t nbk = static_cast<uint8_t>(i - start);

        // Encode super-k-mer: for protein, use 5 bits per AA
        // Calculate bytes needed
        const size_t total_bits = superk_len * 5;
        const size_t nb_bytes = (total_bits + 7) / 8;
        std::vector<uint8_t> encoded(nb_bytes, 0);

        // Encode the super-k-mer sequence
        size_t bit_pos = 0;
        // First k-mer
        for (size_t j = 0; j < kmer_size; ++j)
        {
          const uint8_t code = char_to_code(Alphabet::PROTEIN, kmers[start].at(j));
          // Pack 5 bits into encoded buffer
          for (int b = 4; b >= 0; --b)
          {
            if (code & (1 << b))
              encoded[bit_pos / 8] |= (1 << (7 - (bit_pos % 8)));
            bit_pos++;
          }
        }

        // Remaining bases (last base of each subsequent k-mer)
        for (size_t j = start + 1; j < i; ++j)
        {
          const uint8_t code = char_to_code(Alphabet::PROTEIN, kmers[j].at(kmer_size - 1));
          for (int b = 4; b >= 0; --b)
          {
            if (code & (1 << b))
              encoded[bit_pos / 8] |= (1 << (7 - (bit_pos % 8)));
            bit_pos++;
          }
        }

        // Write to appropriate partition
        writer->insertSuperkmer(encoded.data(), nb_bytes, nbk, part);

        start = i;
      }
    }
  });

  // Cleanup
  delete writer;

  spdlog::info("Protein super-k-mer generation complete for {}", sample_id);
}

} // namespace km
