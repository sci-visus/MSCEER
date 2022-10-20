#ifndef TIMING_H
#define TIMING_H

#include <stdio.h>
#include <vector>
#include <chrono>
#include "gi_basic_types.h"


namespace GInt {

    typedef std::chrono::time_point < std::chrono::high_resolution_clock> TIMEP;
    class ThreadedTimer {
    protected:
        struct Timing {
            int activity;
            TIMEP start_time;
        };
        struct ThreadTimings {
            int thread_id;
            std::vector < Timing > timings;
        };

        ThreadTimings* m_thread_timings;
        int m_num_threads;
        TIMEP m_global_start;
        TIMEP m_global_end;
    public:



        ThreadedTimer(int num_threads) : m_num_threads(num_threads) {
            m_thread_timings = new ThreadTimings[m_num_threads];
        }

        void StartGlobal() {
            m_global_start = std::chrono::high_resolution_clock::now();
        }

        void EndGlobal() {
            m_global_end = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < m_num_threads; i++) this->RecordThreadActivity(i, 0);
        }

        void PrintAll() {
            printf(" took %d ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(m_global_end - m_global_start).count());
        }

        void PrintAllToFile(char* fname, const char** strings) {

            FILE* fout = fopen(fname, "w");
            for (int i = 0; i < m_num_threads; i++) {
                ThreadTimings& tt = m_thread_timings[i];
                for (int a = 0; a < (int)tt.timings.size() - 1; a++) {
                    fprintf(fout, "%d %d %d %s\n", i,
                        std::chrono::duration_cast<std::chrono::milliseconds>(tt.timings[a].start_time - m_global_start).count(),
                        std::chrono::duration_cast<std::chrono::milliseconds>(tt.timings[a + 1].start_time - m_global_start).count(), (char*)(strings[tt.timings[a].activity]));


                }
            }

            fclose(fout);
        }
        void RecordThreadActivity(int thread_num, int activity) {
            Timing t = { activity, std::chrono::high_resolution_clock::now() };
            m_thread_timings[thread_num].timings.push_back(t);
        }

    };

}

#endif
