function plot_TADs(tad_intervals)
    hold on
    for i = 1:length(tad_intervals)-1
        plot([tad_intervals(i)-.5 tad_intervals(i+1)-.5],[tad_intervals(i)-.5 tad_intervals(i)-.5],'k-','linewidth',1)
        plot([tad_intervals(i)-.5 tad_intervals(i+1)-.5],[tad_intervals(i+1)-.5 tad_intervals(i+1)-.5],'k-','linewidth',1)
        plot([tad_intervals(i)-.5 tad_intervals(i)-.5],[tad_intervals(i)-.5 tad_intervals(i+1)-.5],'k-','linewidth',1)
        plot([tad_intervals(i+1)-.5 tad_intervals(i+1)-.5],[tad_intervals(i)-.5 tad_intervals(i+1)-.5],'k-','linewidth',1)
    end
end

