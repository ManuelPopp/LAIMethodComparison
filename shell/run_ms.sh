#!/bin/bash

THRESHOLD=$((250 * 1024 * 1024))  # GB in KB
max_tasks=2

idx_start=1
idx_end=60

# Array to keep track of tasks that need retry
retry_queue=()

# Launch a task function
task() {
    plotid=$1
    log_file="./log/nohup${plotid}ms.log"
    echo "Starting task for plot $plotid"

    Rscript /home/poppman/git/LAIMethodComparison/rsc/multispectral_metrics.R "$plotid"  &> "$log_file" &
    PID=$!

    echo "$PID:$plotid" >> .running_tasks  # store PID and plot
}

# Function to monitor RAM and enforce limit
monitor_ram() {
    while true; do
        if [[ -f .running_tasks ]]; then
            total_rss=0
            youngest_pid=0
            youngest_time=0

            while IFS=: read -r pid plotid; do
                if ps -p "$pid" > /dev/null 2>&1; then
                    rss=$(ps -o rss= -p "$pid" | awk '{sum+=$1} END{print sum}')
                    total_rss=$((total_rss + rss))

                    start_time=$(ps -o etimes= -p "$pid") # elapsed seconds
                    if (( start_time < youngest_time || youngest_time == 0 )); then
                        youngest_time=$start_time
                        youngest_pid=$pid
                    fi
                fi
            done < .running_tasks

            if (( total_rss > THRESHOLD )); then
                echo "RAM usage ${total_rss} KB exceeds threshold ${THRESHOLD} KB. Killing youngest PID $youngest_pid."
                plot_to_retry=$(grep "^$youngest_pid:" .running_tasks | cut -d: -f2)
                kill -9 "$youngest_pid"
                retry_queue+=("$plot_to_retry")
                # remove from .running_tasks
                grep -v "^$youngest_pid:" .running_tasks > .running_tasks.tmp && mv .running_tasks.tmp .running_tasks
            fi
        fi
        sleep 10
        # exit loop if no tasks are running and no retries
        if [[ ! -f .running_tasks ]] && [[ ${#retry_queue[@]} -eq 0 ]]; then
            break
        fi
    done
}

# Clean up any previous temporary file
rm -f .running_tasks

# Main loop
all_plots=($(seq $idx_start $idx_end))

while (( ${#all_plots[@]} > 0 || ${#retry_queue[@]} > 0 )); do
    # enforce max_tasks
    running_count=$(jobs -r | wc -l)
    while (( running_count >= max_tasks )); do
        sleep 5
        running_count=$(jobs -r | wc -l)
    done

    # pick next biome to launch: first retry, else new
    if (( ${#retry_queue[@]} > 0 )); then
        next_plot=${retry_queue[0]}
        retry_queue=("${retry_queue[@]:1}")
    else
        next_plot=${all_plots[0]}
        all_plots=("${all_plots[@]:1}")
    fi

    task "$next_plot"

    sleep 1
done &

# Start RAM monitor in background
monitor_ram &

wait
rm -f .running_tasks
echo "All tasks completed."
