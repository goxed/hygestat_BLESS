AUTOMAKE_OPTIONS = foreign
EXTRA_DIST = \
            src \
            include \
            thirdPartyCodes \
            tools
bin_PROGRAMS = \
             hygestat\
             btt
hygestat_SOURCES=       src/hygestat_windows.cpp\
                        src/barcodes_trimming.cpp\
                        src/bowtie_alignments.cpp\
                        src/chromosome_extraction.cpp\
                        src/compute_cytoband.cpp\
                        src/compute_hypergeometric_tests.cpp\
                        src/compute_hypergeometric_two_samples.cpp\
                        src/compute_sequence_reads.cpp\
                        src/count_hits_bless.cpp\
                        src/count_hits_genomic.cpp\
                        src/data_mapping_fastq.cpp\
                        src/evaluate_mappability.cpp\
                        src/process_fastq.cpp\
                        src/rarefaction.cpp\
                        src/read_config_file.cpp\
                        src/sampling_btt_file.cpp\
                        src/set_output.cpp\
                        src/init_thread.cpp\
                        src/telomere_sequence.cpp\
                        src/two_samples_hygeostat.cpp
btt_SOURCES= src/btt.cpp

install-exec-hook:
		cp -fa tools/*.sh $(DESTDIR)$(bindir)
if FASTQC_HOME_SET
		cp -fa tools/fastqc $(DESTDIR)$(bindir)
endif
if BOWTIE_HOME_SET
		cp -fa tools/bowtie $(DESTDIR)$(bindir)
endif
if SAMSTAT_HOME_SET
		cp -fa tools/samstat $(DESTDIR)$(bindir)
endif
