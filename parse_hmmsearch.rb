#!/usr/bin/env ruby
### License: MIT
### copyright: Yosuke Nishimura (2024)

require 'rake'
Version = "0.1.0"

### [parse_hmmsearch.rb] 
### usage: ruby parse_hmmsearch.rb <hmmsearch output> <protein fasta>
###
### (a) parse only domain-level result and consider domain-level i-Evalue
### (b) return non-overlapping best hit for each protein
### (c) merge multiple domain hits of the same pair of protein and hmm (if not 'non-ovelapped' and 'not-reversed' and 'not-too-distant')

### [ToDo]
### (1) RUBY_VERSION check
### (2) command line option
### (3) consider overlap of domain hit

# {{{ usage
$usage = <<EOS
parse_hmmsearch.rb: parse hmmsearch output and make a hit region table, with merging gap-separated hits into one hit.
Currently, beta version.

[usage]
$ ruby parse_hmmsearch.rb -i <hmmsearch output> -f <protein fasta> -o <output directory> [options]

[input / output]
-i, --hmmsearchout [required]  <hmmsearch output>
-f, --fasta        [required]  <protein fasta>
-o, --outdir       [required]  <output dir>
-H, --hmminfo                  <tsv of hmms>

[optional output]
--create-evalue-table [off]   -- create evalue table of all combinations of proteins and hmms

[paramters]
-ge, --gene-evalue  [1e-2]        <gene-level i-Evalue threshold>
-e, --evalue  [10.0]              <domain-level i-Evalue threshold>
--min-hmm-len [0](int)            <minimum length of hmmlen of domain hit>
--min-hmm-cov [0](float, 0-1)     <minimum fraction of hmmlen of domain hit in hmm length>
--min-ali-len [0](int)            <minimum length of alilen of domain hit>
--min-ali-cov [0](float, 0-1)     <minimum fraction of alilen of domain hit in protein length>
--max-ali-ovp-frc [0.2](float, 0-1)   <maximum fraction of overlaping of domain hits to be reported in best-hmm.tsv>
--max-ali-ovp-len [100000](int)       <maximum length of overlaping of domain hits to be reported in best-hmm.tsv>
--max-hmm-ovp-frc [0.2](float, 0-1)   <maximum fraction of overlaping of domain hits to be reported in best-hmm.tsv>
--max-hmm-ovp-len [100000](int)       <maximum length of overlaping of domain hits to be reported in best-hmm.tsv>

[others]
-v, --version
-h, --help
EOS
# }}}

# {{{ threshold / option default
e_thre    = "10.0"
e_gene    = "1e-2"
minhmmlen = "0"
minhmmcov = "0.0"
minalilen = "0"
minalicov = "0.0"
maxaliovpfrc = "0.2"
maxaliovplen = "100000"
maxhmmovpfrc = "0.2"
maxhmmovplen = "100000"
fa = nil ### fasta
fh = nil ### hmmout
odir = nil
hi = nil
is_fa_gz = false
is_fh_gz = false
e_table  = false

### parse command line
i = 0
while i < ARGV.size 
  e = ARGV[i]

  case e
  when "-h", "--help"
    puts $usage
    exit(0)
  when "-v", "--version"
    puts Version
    exit(0)
  when "-i", "--hmmsearchout"
    fh      = ARGV[i+1]; i += 1
  when "-f", "--fasta"
    fa      = ARGV[i+1]; i += 1
  when "-o", "--outdir"
    odir    = ARGV[i+1]; i += 1
  when "-H", "--hmminfo"
    hi = ARGV[i+1]; i += 1
  when "-e", "--evalue"
    e_thre  = ARGV[i+1]; i += 1
  when "-ge", "--gene-evalue"
    e_gene  = ARGV[i+1]; i += 1
  when "--create-evalue-table"
    e_table = true
  when "--min-hmm-len"
    minhmmlen = ARGV[i+1]; i += 1
  when "--min-hmm-cov"
    minhmmcov = ARGV[i+1]; i += 1
  when "--min-ali-len"
    minalilen = ARGV[i+1]; i += 1
  when "--min-ali-cov"
    minalicov = ARGV[i+1]; i += 1
  when "--max-ali-ovp-frc"
    maxaliovpfrc = ARGV[i+1]; i += 1
  when "--max-ali-ovp-len"
    maxaliovplen = ARGV[i+1]; i += 1
  when "--max-hmm-ovp-frc"
    maxhmmovpfrc = ARGV[i+1]; i += 1
  when "--max-hmm-ovp-len"
    maxhmmovplen = ARGV[i+1]; i += 1
  end
  i += 1
end
# }}}

# {{{ input validation
raise("<hmmsearch output> file is not specified.") unless fh
raise("<protein fasta> file is not specified.")    unless fa
raise("<output directory> is not specified.")      unless odir
raise("<hmmsearch output> file does not exist.")   unless File.exist?(fh)
raise("<protein fasta> file does not exist.")      unless File.exist?(fa)
open(fh){ |fr|
  if fh =~ /\.gz$/ or fh =~ /\.gzip/
    is_fh_gz = true
    require 'zlib'

    Zlib::GzipReader.open(fh) { |fr|
      l = fr.gets
      raise("<hmmsearch output> file has wrong format.\nThe first line should be: '# hmmsearch :: search profile(s) against a sequence database'") if l !~ /^# hmmsearch :: /
    }
  else
    l = fr.gets
    raise("<hmmsearch output> file has wrong format.\nThe first line should be: '# hmmsearch :: search profile(s) against a sequence database'") if l !~ /^# hmmsearch :: /
  end
}
open(fa){ |fr|
  if fa =~ /\.gz$/ or fa =~ /\.gzip/
    is_fa_gz = true
    require 'zlib'

    Zlib::GzipReader.open(fa) { |fr|
      l = fr.gets
      raise("<protein fasta> file has wrong format.\nThe first line should begin '>'") if l[0] != ">"
    }
  else
    l = fr.gets
    raise("<protein fasta> file has wrong format.\nThe first line should begin '>'") if l[0] != ">"
  end
}
# }}}

### float checker (https://stackoverflow.com/questions/1034418/determine-if-a-string-is-a-valid-float-value)
class String
  def is_float?
    # The double negation turns this into an actual boolean true - if you're 
    # okay with "truthy" values (like 0.0), you can remove it.
    !!Float(self) rescue false
  end
  def is_integer?
    self.to_i.to_s == self
  end
end
[e_thre, minhmmcov, minalicov, maxaliovpfrc, maxhmmovpfrc].zip(
  %w|evalue min-hmm-cov min-ali-cov max-ali-ovp-frc max-hmm-ovp-frc|){ |var, str|
  raise("'--#{str} #{var}' is not valid. Must be float") unless var.is_float?
}
[minhmmlen, minalilen, maxaliovplen, maxhmmovplen].zip(
  %w|min-hmm-len min-ali-len max-ali-ovp-len max-hmm-ovp-len|){ |var, str|
  raise("'--#{str} #{var}' is not valid. Must be integer") unless var.is_integer?
}

E_thre    = e_thre.to_f
E_gene    = e_gene.to_f
E_table   = e_table
MinHmmLen = minhmmlen.to_i
MinAliLen = minalilen.to_i
MinHmmCov = minhmmcov.to_f
MinAliCov = minalicov.to_f
MaxAliOvpFrc = maxaliovpfrc.to_f
MaxAliOvpLen = maxaliovplen.to_i
MaxHmmOvpFrc = maxhmmovpfrc.to_f
MaxHmmOvpLen = maxhmmovplen.to_i

### make output directory
$stderr.puts "output directory: #{odir} exists. Overwrite it." if File.directory?(odir)
require 'fileutils'
FileUtils.mkdir_p(odir)


# {{{ [0] range function
class Range
  include Comparable

  def <=>(other)
    self.min <=> other.min
  end
  def overlap?(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return false
    else
      return true
    end
  end
  def &(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return nil
    else
      return ([s.min, t.min].max)..([s.max, t.max].min)
    end
  end
  def |(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return [s, t]
    else
      return [([s.min, t.min].min)..([s.max, t.max].max)]
    end
  end
  def include?(t)
    s = self
    o = s & t  ## overlap of s and t
    return false unless o
    if o.min == t.min and o.max == t.max ## s includes t
      return true
    else
      return false
    end
  end
end
def merge_ranges(ranges) ## ranges = [3..300, 320..500, 504..732, ...]
  return ranges if ranges.size < 2

  rs  = ranges.sort
  (-rs.size..-2).each{ |j| ## index from right
    merged = rs[j] | rs[j+1]
    case merged.size
    when 1 ## overlap detected
      rs[j] = merged[0]
      rs.delete_at(j+1)
    when 2 ## overlap not detected
      ## do nothing
    else raise
    end
  }
  return rs
end
# }}} range function


# {{{ [1] parse gene ids from hmmsearchout
def parse_hmmsearch_1(fh, is_fh_gz, _gids)
  ### scan only gene id
  fr = is_fh_gz ? Zlib::GzipReader.open(fh) : open(fh)

  while l = fr.gets
    if l =~ /^Scores for complete sequences/
      flag = "parse_full"
    elsif l =~ /--- inclusion threshold ---/
      flag = "above_inclusion_threshold"
    elsif l =~ /^$/ ## reset parse flag if a blank line comes
      flag = ""
    elsif flag == "parse_full" and l =~ /^\s+\d/
      full_e, score, gid = l.strip.split(/\s+/).values_at(0, 1, 8) ## full evalue

      ### full evalue cutoff
      next if full_e.to_f >= E_gene

      _gids[gid] = 1
    end

  end
end
_gids = {}
parse_hmmsearch_1(fh, is_fh_gz, _gids)
$stderr.puts "hit gene (n=#{_gids.size} entries) parsed from hmmsearchout. [#{Time.now}]"
# }}}


# {{{ [2] parse fasta, gene ids (ordered by a given fasta file)
gids      = {}
gid2info  = {}
gid2seq   = {}
gid2len   = {}

def parse_fasta(fa, is_fa_gz, _gids, gids, gid2info, gid2seq, gid2len)
  fr = is_fa_gz ? Zlib::GzipReader.open(fa) : open(fa)

  seq       = ""
  gid       = ""
  while l = fr.gets
    if l[0] == ">"
      ### store seq
      if _gids[gid]
        gid2seq[gid] = seq
        gid2len[gid] = seq.size
      end

      ### init
      gid, *info = l[1..-1].split(/\s/)
      seq = ""

      ### store gid and header
      if _gids[gid]
        gids[gid] = 1
        gid2info[gid] = info*" "
      end
    else
      seq += l.strip if _gids[gid]
    end
  end

  ### store seq
  if _gids[gid]
    gid2seq[gid] = seq
    gid2len[gid] = seq.size
  end

  fr.close
end
parse_fasta(fa, is_fa_gz, _gids, gids, gid2info, gid2seq, gid2len)
$stderr.puts "protein fasta (n=#{gids.size} entries) parsed. [#{Time.now}]"
# }}}


### [3] parse hmmsearch output
hmms       = {}
evalues    = {} ## store domain i-Evalue
all_info   = {} ## store all column
ranks      = {} ## store rank based on evalue
best_idxs  = {} ## store best index based on evalue
hmm2acc    = {}
hmm2desc   = {}
hmm2len    = {}

fr = is_fh_gz ? Zlib::GzipReader.open(fh) : open(fh)
# Query:       1-cysPrx_C  [M=40]
# Accession:   PF10417.9
# Description: C-terminal domain of 1-Cys peroxiredoxin
# Scores for complete sequences (score includes all domains):
#    --- full sequence ---   --- best 1 domain ---    -#dom-
#     E-value  score  bias    E-value  score  bias    exp  N  Sequence  Description
#     ------- ------ -----    ------- ------ -----   ---- --  --------  -----------
#     4.6e-12   51.2   0.4    6.3e-12   50.7   0.4    1.2  1  gene1
#       0.004   22.6   0.5      0.016   20.6   0.0    2.4  1  geneX
#   ------ inclusion threshold ------
#       0.011   21.2   0.0       0.03   19.7   0.0    1.8  1  geneY

# >> geneX
#    #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
#  ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
#    1 ?    6.6   0.0    0.0083   2.1e+02      90     135 ..      16      57 ..       5      74 .. 0.79
#    2 ?    7.7   0.0     0.004     1e+02      89     134 ..     123     164 ..      94     173 .. 0.80
#
flag  = ""
gid   = ""
hmm   = ""
n     = 0
inc_gids = {} ### gene ids above inclusion threshold

while l = fr.gets
  if l =~ /^Query:\s+(\S+)\s+\[M=(\d+)\]/  # detect new entry
    hmm, hlen = $1, $2
    hmms[hmm] = 1
    hmm2len[hmm] = hlen.to_i

  elsif l =~ /^Accession:\s+(\S+)/
    acc    = $1.strip
    hmm2acc[hmm] = acc

  elsif l =~ /^Description:\s+(.+)/
    desc   = $1.strip
    hmm2desc[hmm] = desc

  elsif l =~ /^Scores for complete sequences/
    flag   = "parse_full"

  elsif l =~ /--- inclusion threshold ---/
    flag   = "above_inclusion_threshold"

  elsif l =~ /^$/ ## reset parse flag if a blank line comes
    flag   = ""

  elsif flag == "parse_full" and l =~ /^\s+\d/
    full_e, score, gid = l.strip.split(/\s+/).values_at(0, 1, 8) ## full evalue

    ### full evalue cutoff
    next if full_e.to_f >= E_gene

    inc_gids[gid]      ||= {}
    inc_gids[gid][hmm]   = [full_e.to_f, score.to_f]

  elsif flag == "above_inclusion_threshold" and l =~ /^\s+\d/
    # full_e, gid = l.strip.split(/\s+/).values_at(0, 8) ## full evalue

  elsif l =~ /^>>\s+(\S+)/ ## reset parse flag if a blank line comes
    gid  = $1

    next unless inc_gids[gid]      ### [!!!] skip if below_inclusion_threshold
    next unless inc_gids[gid][hmm] ### [!!!] skip if below_inclusion_threshold

    gE = inc_gids[gid][hmm][0] ### gene-level evalue
    next if gE >= E_gene ### [!!!] skip if above gene-level evalue threshold

    flag = "parse_each"

    ### abort if gid is not found in fasta.
    # raise("#{fh}: #{gid} is reported, but not found in #{fa}") unless gids[gid]

  elsif flag == "parse_each" and l =~ /^\s+\d/
    a = l.strip.split(/\s+/)
    evalue, hmm_fm, hmm_to, ali_fm, ali_to = a.values_at(5, 6, 7, 9, 10)
    hmmlen = hmm_to.to_i - hmm_fm.to_i + 1
    alilen = ali_to.to_i - ali_fm.to_i + 1

    ### check various threshold 
    next if evalue.to_f >= E_thre
    next if hmmlen < MinHmmLen
    next if alilen < MinAliLen
    next if hmmlen.to_f / hmm2len[hmm] < MinHmmCov
    next if alilen.to_f / gid2len[gid] < MinAliCov

    evalues[gid]       ||= {}
    evalues[gid][hmm]  ||= []
    evalues[gid][hmm]  << evalue.to_f

    all_info[gid]      ||= {}
    all_info[gid][hmm] ||= []

    info  = a.values_at(2..7, 9..10, 12..13, 15) ### score, bias, c-Evalue, i-Evalue, hmmfrom, hmmto, alifrom, alito, envfrom, envto, acc
    info += inc_gids[gid][hmm] ### [full evalue, full score]
    all_info[gid][hmm] << info
  end

  n += 1
  $stderr.puts "parse hmmsearch out: lines #{n} [#{Time.now}]" if n % 5_000_000 == 0
end

fr.close
$stderr.puts "hmmsearch out parsed. [#{Time.now}]"


# {{{ E_table
if E_table
  ### fill ranks and best_idxs
  gids.each_key.with_index{ |gid, n|
    $stderr.puts "best evalue detection: #{n}th: #{gid} [#{Time.now}]" if n % 50000 == 0 and n > 0

    next unless evalues[gid]

    hmms.each_key{ |hmm|
      e = evalues[gid][hmm] ### array of evalues
      next unless e

      r = e.map.with_index.sort.map(&:last) ## https://stackoverflow.com/questions/14446181/originals-indexes-of-sorted-elements-in-ruby
      ranks[gid]    ||= {}
      ranks[gid][hmm] = r

      bidx = nil
      r.each.with_index{ |i, idx| bidx = idx if i == 0 }

      best_idxs[gid]    ||= {}
      best_idxs[gid][hmm] = bidx
    }
  }
  $stderr.puts "best_idxs parsed. [#{Time.now}]"
  # require 'pp' ; pp best_idxs


  ### evalues.tsv
  open("#{odir}/evalues.tsv", "w"){ |fw|
    header  = ["seq"]
    header += hmms.keys.map{ |hmm| hmm2acc[hmm] || hmm }
    fw.puts header*"\t"

    gids.each_key{ |gid|
      o = [gid]

      if evalues[gid]
        hmms.each_key{ |hmm|
          if evalues[gid][hmm]
            bidx = best_idxs[gid][hmm]
            o << evalues[gid][hmm][bidx]
          else
            o << ""
          end
        }
      else
        hmms.size.times do o << "" end
      end

      fw.puts o*"\t"
    }
  }
  $stderr.puts "evalues.tsv finished. [#{Time.now}]"
  # require 'pp' ; pp evalues
end
# }}}


### all  (all combination and all domain hit)
fwf  = open("#{odir}/all-hit.tsv", "w")
### best (all combination and best domain hit considering overlap)
fwb  = open("#{odir}/best-hit.tsv", "w")
### best fasta
fwbW = open("#{odir}/best-hit.whole.fa", "w")
fwbF = open("#{odir}/best-hit.fa", "w")

fws = [fwf, fwb, fwbW, fwbF] ### writers

# bhmm = {} ## bhmm[hmm] = [gid, evalue] ### best gid for each hmm
# bgid = {} ## bgid[gid] = [hmm, evalue] ### best hmm for each gid

### header
# header = %w|protein length(aa) protein_info hmm_name hmm_acc hmm_desc hmm_len score bias c-Evalue i-Evalue hmm.fm hmm.to ali.fm ali.to env.fm env.to acc|
header = %w|protein length(aa) protein_info hmm_name hmm_acc hmm_desc hmm_len score bias c-Evalue i-Evalue hmm.fm hmm.to ali.fm ali.to env.fm env.to acc full-Evalue full-score|
fwf.puts header*"\t"
fwb.puts (header + %w|link region_name|)*"\t"

gids.each_key{ |gid|
  next unless evalues[gid]

  ginfo = gid2info[gid]
  glen  = gid2len[gid]
  # gseq  = gid2seq[gid]

  infos  = [] ### all  info (array of info)
  alis   = [] ### [31..201, 300..403, ..] ### store non-overllapping hit range of ali 
  binfos = [] ### best info (array of info)

  ### store all hit for this gid into infos
  hmms.each_key{ |hmm|
    next unless evalues[gid][hmm]
    all_info[gid][hmm].each{ |a|
      ary = [hmm] + a
      infos << ary
    }
  }

  ### sort all hit for this gid
  infos = infos.sort_by{ |info| info[4].to_f } ### sorted by i-Evalue

  ### store best hit into binfos
  infos.each{ |info|
    ### [hmm, score, bias, c-Evalue, i-Evalue, hmmfrom, hmmto, alifrom, alito, envfrom, envto, acc]
    ###    0,     1,    2,,       3,        4,       5,     6,       7,     8,       9,    10,  11]
    ali_fm = info[7].to_i
    ali_to = info[8].to_i
    ali = ali_fm..ali_to
    ali_len = ali.max - ali.min + 1

    flag = 0 ### overlap flag
    alis.each{ |_ali|
      _ali_len = _ali.max - _ali.min + 1

      intersect = _ali & ali
      if intersect ### overlapping
        ovp_len = intersect.max - intersect.min + 1
        ovp_frc = ovp_len.to_f / [ali_len, _ali_len].min

        ### overlapped
        if ovp_len > MaxAliOvpLen or ovp_frc > MaxAliOvpFrc
          flag = 1
        end
      end
    }

    ### add ali and best_info
    if flag == 0
      alis   << ali
      binfos << info
    end
  }

  ### make full output
  infos = infos.sort_by{ |info| [info[7].to_i, info[4].to_f] } ### ali.fm, i-Evalue
  infos.each{ |info|
    hmm = info[0]
    pref = [gid, glen, ginfo, hmm, hmm2acc[hmm], hmm2desc[hmm], hmm2len[hmm]]

    o = pref + info[1..-1]

    fwf.puts o*"\t"
  }

  ### make best output
  binfos = binfos.sort_by{ |info| [info[7].to_i, info[4].to_f] } ### ali.fm, i-Evalue

  ## find links (the same pair of protein and hmm (if not 'non-ovelapped' and 'not-reversed'))

  # (link1) Arc_Ant_N02_N0000523_26 1049  TIGR03479 DMSO_red_II_alp: DMSO reductase family type II enzyme, molybdopterin subunit  913 28.2  0 7.60E-11  2.50E-05  104 150 154 200 148 203 0.96
  # (link1) Arc_Ant_N02_N0000523_26 1049  TIGR03479 DMSO_red_II_alp: DMSO reductase family type II enzyme, molybdopterin subunit  913 37.3  0 1.30E-13  4.50E-08  252 463 332 561 273 565 0.78
  #         Arc_Ant_N02_N0000523_26 1049  TIGR00509 bisC_fam: molybdopterin guanine dinucleotide-containing S/N-oxide reductases  772 29.1  0 4.50E-11  1.50E-05  601 706 876 982 861 991 0.82
  # (link2) Arc_Ant_N02_N0000547_26 354   TIGR03408 urea_trans_UrtC: urea ABC transporter, permease protein UrtC  313 28.1  4.2   1.60E-09  8.50E-05   8  53  39   84 32   90 0.93
  # (link2) Arc_Ant_N02_N0000547_26 354   TIGR03408 urea_trans_UrtC: urea ABC transporter, permease protein UrtC  313 93.9  19.3  1.50E-29  7.90E-25  88  312 95  346 88  347 0.81

  link_cands = {}

  binfos.each.with_index{ |info, idx|
    hmm = info[0]
    link_cands[hmm] ||= []
    link_cands[hmm] << (info + [idx])
  }

  idx2link  = {} ### h[0] = 1, h[1] = 1, h[2] = 1, h[3] = 2, h[4] = 2 ### 0-based --> 1-based
  idx2range = {} ### h[0] = 10..230, h[1] = 10..230, h[2] = 10..230, h[3] = 450..540, h[4] = 450..540 ### merged ali range
  link_cands.each_key{ |hmm|
    pairs = [] ### pairs[0] = 0..1, pairs[1] = 1..2, pairs[2] = 3..4
    links = [] ### links[0] = 0..2, links[1] = 3..4

    cands = link_cands[hmm]
    n = cands.size
    next if n <= 1

    ### store pairs
    (0..n-2).each{ |i|
      i_fm = cands[i][5].to_i ### hmm from
      i_to = cands[i][6].to_i ### hmm to
      i_range = i_fm..i_to
      i_len   = i_to - i_fm + 1

      j = i+1
      j_fm = cands[j][5].to_i ### hmm from
      j_to = cands[j][6].to_i ### hmm to
      j_range = j_fm..j_to
      j_len   = j_to - j_fm + 1

      ### [!!!] SKIP if hmm from is reversed 
      next if i_fm > j_fm

      ### non-overlapped
      intersect = i_range & j_range
      flag = 0 ### flag = 1 if overlapped
      if intersect
        ovp_len = intersect.max - intersect.min + 1
        ovp_frc = ovp_len.to_f / [i_len, j_len].min
        if ovp_len > MaxAliOvpLen or ovp_frc > MaxAliOvpFrc
          flag = 1
        end
      end

      if flag == 0 ### non-overlapped
        pairs << (i..j)
      end
    }

    ### store links
    links = merge_ranges(pairs)
    links.each.with_index{ |link, idx|
      elems = link.to_a ## 0..2 --> [0, 1, 2]
      fm = elems.map{ |i| cands[i][7].to_i }.min ### min of ali.fm
      to = elems.map{ |i| cands[i][8].to_i }.max ### max of ali.to

      elems.each{ |i|
        ### i: idx for link_cands (each hmm)
        ### j: idx for binfos
        j = cands[i][-1]

        idx2link[j]  = idx + 1
        idx2range[j] = fm..to
      }
    }
  }

  ### write output
  binfos.each.with_index{ |info, idx|
    hmm = info[0]
    pref = [gid, glen, ginfo, hmm, hmm2acc[hmm], hmm2desc[hmm], hmm2len[hmm]]

    ali_fm = info[7].to_i
    ali_to = info[8].to_i

    if link = idx2link[idx]
      fm = idx2range[idx].min
      to = idx2range[idx].max
      name   = "#{gid}_fm#{fm}_to#{to}"
      o = pref + info[1..-1] + ["link#{link}_#{hmm}", name]
    else
      name   = "#{gid}_fm#{ali_fm}_to#{ali_to}"
      o = pref + info[1..-1] + ["", name]
    end

    fwb.puts o*"\t"

    name =~ /_fm(\d+)_to(\d+)$/
    fm = $1.to_i - 1
    to = $2.to_i - 1

    labF = "#{name} #{gid2info[gid]}".strip
    labW = "#{gid} #{gid2info[gid]}".strip
    fwbF.puts ">#{labF}\n#{gid2seq[gid][fm..to]}"
    fwbW.puts ">#{labW}\n#{gid2seq[gid]}"
  }
}

### file close
fws.each{ |fw| fw.close }
$stderr.puts "all process finished. [#{Time.now}]"
