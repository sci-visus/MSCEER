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
                for (int a = 0; a < tt.timings.size() - 1; a++) {
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

	class SimpleTimer {
	public:
		enum printmode {
			TIMER_PRINT_START_END,
			TIMER_PRINT_AT_TASK_END,
			TIMER_PRINT_ALL_AT_END,
			TIMER_SILENT
		};
	protected:
		typedef std::chrono::steady_clock::time_point TPOINT;
		TPOINT m_global_start_time;
		bool m_has_global_start_time;
		struct Timing {
			std::string taskname;
			TPOINT start;
			TPOINT end;
		};

		std::vector<Timing> m_timings_sequence;

	protected:
		printmode m_mode;
	public:
		SimpleTimer(printmode m) : m_mode(m), m_has_global_start_time(false) {
			m_timings_sequence.reserve(100);
		}
		
		void StartGlobal() {
			m_global_start_time = std::chrono::steady_clock::now();
			m_has_global_start_time = true;
			StartTask("Timer Global Start");
		}

		// returns task id
		int StartTask(std::string s) {
			if (!m_has_global_start_time) StartGlobal(); // start the global start time if user did not do it
			TPOINT task_start_time = std::chrono::steady_clock::now();
			m_timings_sequence.push_back({ s, task_start_time, task_start_time });
			
			if (m_mode == TIMER_PRINT_START_END) {
				printf("TIMER: Starting Task %s at %d ms\n", s.c_str(),
					std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - m_global_start_time).count());
			}
			return m_timings_sequence.size() - 1;
		}

		// use task id to record the end of a task - will print depending on mode
		void EndTask(int task_id) {
			TPOINT now_time = std::chrono::steady_clock::now();
			Timing& t = m_timings_sequence[task_id];
			t.end = now_time;

			if (m_mode == TIMER_PRINT_START_END ||
				m_mode == TIMER_PRINT_AT_TASK_END) {
				// format: [global activity name] [task] [start] [end] [dration]

				printf("TIMER: TASK: %s START: %d END: %d DURATION: %d\n", t.taskname.c_str(),
					std::chrono::duration_cast<std::chrono::milliseconds>(t.start - m_global_start_time).count(),
					std::chrono::duration_cast<std::chrono::milliseconds>(t.end - m_global_start_time).count(),
					std::chrono::duration_cast<std::chrono::milliseconds>(t.end - t.start).count());
			}
		}

		void EndGlobal() {
			EndTask(0); // global task was the first one we added

			if (m_mode == TIMER_PRINT_ALL_AT_END) {
				for (auto t : m_timings_sequence) {
					printf("TIMER: TASK: %s START: %d END: %d DURATION: %d\n", t.taskname.c_str(),
						std::chrono::duration_cast<std::chrono::milliseconds>(t.start - m_global_start_time).count(),
						std::chrono::duration_cast<std::chrono::milliseconds>(t.end - m_global_start_time).count(),
						std::chrono::duration_cast<std::chrono::milliseconds>(t.end - t.start).count());
				}
			}

		}

		void WriteTimingsToFile(std::string fname) {
			FILE* fout = fopen(fname.c_str(), "w");
			for (auto t : m_timings_sequence) {
				fprintf(fout, "TIMER: TASK: %s START: %d END: %d DURATION: %d\n", t.taskname.c_str(),
					std::chrono::duration_cast<std::chrono::milliseconds>(t.start - m_global_start_time).count(),
					std::chrono::duration_cast<std::chrono::milliseconds>(t.end - m_global_start_time).count(),
					std::chrono::duration_cast<std::chrono::milliseconds>(t.end - t.start).count());
			}
			fclose(fout);
		}
	};

}

#endif
