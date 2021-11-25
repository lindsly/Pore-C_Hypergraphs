function hfchrom = chromosomeplot(cytobandstruct, chromnum, varargin)
% CHROMOSOMEPLOT Plots chromosome ideograms with the G-banding pattern from Giemsa staining.
%
% CHROMOSOMEPLOT(S) plots the ideogram of all chromosomes (summary view)
% using the G-Banding pattern information from S. S can be a structure
% containing the fields ChromLabels, BandStartBPs, BandEndBPs, BandLabels,
% and GieStains, such as returned by CYTOBANDREAD function. S can also be
% the file name of an NCBI ideogram text file or a UCSC Genome Browser
% cytoband text file. The Giemsa staining produces G-bands indicating
% different densities (gneg, gpos25, gpos50, gpos75 and gpos100),
% centromeres (acen), repeats (stalk), and constitutive heterochromatins
% (gvar). 
% 
% CHROMOSOMEPLOT(S,CHROMNUM) plots the ideogram of a single chromosome
% specified by CHROMNUM. CHROMNUM must be one of the chromosome numbers in
% the genome of structure S. It can be a numeric number or a string. For
% example, 15 for chromosome 15, and 'X', for chromosome X. If CHROMNUM is
% set to 0, then all the chromosome ideograms will be displayed.
% 
% CHROMOSOMEPLOT(..., 'ORIENTATION', ORIEN) specifies the orientation of
% the ideogram when plotting a single chromosome. The choices are
% 'Vertical' or 1 and 'Horizontal' or 2. The default is 'Vertical'. This
% option is ignored in summary view, which displays all the chromosomes
% vertically.
% 
% CHROMOSOMEPLOT(..., 'SHOWBANDLABEL', TF) displays band labels when
% plotting a single chromosome.  Default is TRUE. Band labels will not be
% displayed in summary view.
% 
% CHROMOSOMEPLOT(..., 'ADDTOPLOT', AXESHANDLE) plots a single chromosome
% ideogram specified by CHROMNUM to a figure axis AXESHANDLE. This option
% is ignored in summary view.
% 
% CHROMOSOMEPLOT(..., 'UNIT', U) specifies the units (base pair, kilo base
% pair, and mega base pair) used to display the start and end genome
% positions in the data tip of the ideogram. Choices are 1(BP), 2(KB), or
% 3(MB).  
% 
% CHROMOSOMEPLOT(..., 'CNV', CNVSTRUCT) displays copy number variance (CNV)
% data aligned to the chromosomes in the ideogram. Gains are shown in green
% to the right or above the ideogram. Losses are shown in red to the left
% or below the ideogram. The CNVSTRUCT is a structure array containing
% these fields:
%       Chromosome - Chromosome numbers 
%       CNVType - Loss (1), Gain (2)
%       Start - Start genomic position of the CNV (in BP units)
%       End - End genomic position of the CNV (in BP units)  
% 
% Note: The Chromosome field can be a numeric or character vector of
% chromosome numbers. 
% 
%   Examples:
%
%       % Read the cytogenetic banding information for Homo sapiens and 
%       % plot the chromosome ideograms of Homo sapiens
%       hs_cytobands = cytobandread('hs_cytoBand.txt')   
%       chromosomeplot(hs_cytobands);
%       title('Human Karyogram')
%
%       % Plot chromosome 15 of Homo sapiens horizontally and display 
%       % the genomic positions in the datatip in KB units 
%       chromosomeplot(hs_cytobands, 15, 'Orientation', 2, 'Unit', 2);
% 
%       % Add chromosome 10 ideogram to a plot of array CGH data of
%       % chromosome 10 from sample #3 in the Coriell cell line study. 
%       % The genomic position in this study is in kilo base pair unit.  
%       load coriell_baccgh
%       S = cghcbs(coriell_data,'sampleind',3,'chromosome',10,'showplot',10);
%       chromosomeplot('hs_cytoBand.txt', 10, 'addtoplot', gca, 'unit', 2)
% 
%       % Plot an CNV data aligned to chromosomes in human ideogram
%       cnvStruct = struct('Chromosome', char({'10', 'X'}),...
%                          'CNVType', [2 1],...
%                          'Start', [66905000 25416000],...
%                          'End',   [110412000 53357000]);
%       chromosomeplot('hs_cytoBand.txt', 'cnv', cnvStruct);
%
%   See also AFFYSNPCNVDEMO, BACACGHDEMO, CGHCBS, CYTOBANDREAD.

% CHROMOSOMEPLOT(..., 'FIGURE', FIGUERHANDLE) plots chromosome(s) to
% specified figure FIGUREHANDLE.

%   Copyright 2007-2011 The MathWorks, Inc.



% Input checking
if nargin > 0
    cytobandstruct = convertStringsToChars(cytobandstruct);
end

if nargin > 1
    chromnum = convertStringsToChars(chromnum);
end

if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

bioinfochecknargin(nargin,1,mfilename);

if ~isstruct(cytobandstruct)
    try
        cytobandstruct = cytobandread(cytobandstruct);
    catch theException
        error(message('bioinfo:chromosomeplot:UnableReadInputFile', theException.message));
    end
elseif ~checkCytoBandStruct(cytobandstruct)
   error(message('bioinfo:chromosomeplot:BadFieldNames'));
end

if  nargin == 1 
    chromnum = 0;  % will plot all the chromosomes
end

if nargin >= 2
    if isnumeric(chromnum) && isscalar(chromnum)
        if ~isreal(chromnum) || chromnum < 0
            error(message('bioinfo:chromosomeplot:ChromNumNotANumber'));
        end
        chromnum = fix(chromnum);
    elseif ischar(chromnum)
        tmp = str2double(chromnum); %ok
        if ~isnan(tmp)
            chromnum = tmp;
        end
    else
        error(message('bioinfo:chromosomeplot:ChromNumNotCorrectType'));
    end
end

% Check the number of chromosomes in cytobandstruct
chromLabels = unique(cytobandstruct.ChromLabels);
tmparg = chromnum;
chromnum = validateChromNumber(chromLabels, chromnum);
if chromnum < 0 && rem(nargin,2) ~= 1
    error(message('bioinfo:chromosomeplot:ChromNumNotValid'));
end

% Initialization
orientation = 1;
hf_old = [];
ha_old = [];
addtoplot = false;
showbandlabels = true;
cnvstruct = [];
unit = 1; % 1 - 1B, 2-KB, 3-MB

if nargin > 2
    if chromnum < 0;
        if rem(nargin,2) == 0
            error(message('bioinfo:chromosomeplot:IncorrectNumberOfArguments', mfilename));
        end
        pvnum = 1; % Only struct input
        cpargs = {tmparg, varargin{:}}; %#ok<*CCAT>
        chromnum = 0;
    else
        if rem(nargin,2) == 1
            error(message('bioinfo:chromosomeplot:IncorrectNumberOfArguments', mfilename));
        end
        pvnum = 2; % struct+chromnum input
        cpargs = varargin;
    end
    
    okargs = {'orientation', 'addtoplot', 'cnv', 'showbandlabel','unit','figure'};
    for j=1:2:nargin - pvnum
        pname = cpargs{j};
        pval = cpargs{j+1};
        k = find(strncmpi(pname, okargs, numel(pname)));
        if isempty(k)
            error(message('bioinfo:chromosomeplot:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:chromosomeplot:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1
                    if ischar(pval) 
                        if strncmpi(pval, 'vertical', 1)
                            orientation = 1;
                        elseif strncmpi(pval, 'horizontal', 1)
                            orientation = 2;
                        else
                            warning(message('bioinfo:chromosomeplot:InvalidOrientation', pval));
                        end
                    elseif isscalar(pval)
                        if pval ~= 1 && pval ~= 2
                            warning(message('bioinfo:chromosomeplot:InvalidNumericOrientation', pval));
                        else
                            orientation = pval;
                        end
                    end
                case 2
                    if ishandle(pval) && strcmpi(get(pval, 'type'), 'axes')
                        ha_old = pval;
                        addtoplot = true;
                    else
                        error(message('bioinfo:chromosomeplot:InputNotAnAxesHandle'));
                    end
                case 3
                    cnvstruct = pval;
                    if ~isstruct(cnvstruct) || ~all(size(cnvstruct)==[1,1])|| ~checkCNVStruct(cnvstruct) 
                        error(message('bioinfo:chromosomeplot:CNVInfoNotValid'));
                    end
                case 4  % s flag
                    showbandlabels = bioinfoprivate.opttf(pval,okargs{k},mfilename);              
                case 5 % unit
                    if isnumeric(pval) && isscalar(pval) && pval >=1 && pval <=3
                        unit = fix(pval);
                    else
                         error(message('bioinfo:chromosomeplot:IncorrectUnit'));
                    end
                    
                case 6 % figure
                    if ishandle(pval) && strcmpi(get(pval, 'type'), 'figure')
                        hf_old = pval;
                    else
                        error(message('bioinfo:chromosomeplot:InputNotAFigureHandle'));
                    end
            end
        end
    end
end

%  Only show band labels in sigle chromosome mode. Not in summary mode
if chromnum == 0
    showbandlabels = false;
    
    if orientation == 2
        orientation = 1;
        warning(message('bioinfo:chromosomeplot:OnlyVerticalSummaryView'));
    end
end

%--------Done input variable checking--------------------------%
chromaxes_props = { 'Tag', 'chromplotaxes',...
                    'Box', 'off',... 
                    'Xtick',[],'Ytick',[],...
                    'Xlim', [0,1],...
                    'Ylim', [0,1],...
                    'Color', [1, 1, 1]};%,...

% Plot all the chromosomes in one view
if addtoplot
    if chromnum == 0
        error(message('bioinfo:chromosomeplot:ChromNumNotSpecified'));
    end
    
    hf = get(ha_old, 'Parent');
    ch = get(hf, 'Children');
    axs = findobj(ch, 'Type', 'axes');

    if numel(axs) > 1
        error(message('bioinfo:chromosomeplot:TooManyAxesinFigure'));
    end

    h = 0.15; % Height for plot the chromosome
    voffset = 0.05; % vertical offset from the bottom
    % Reposition the original axes and add a new axes
    rect_o = get(ha_old, 'Position');
    set(ha_old, 'Position', [rect_o(1) h+voffset rect_o(3) rect_o(4)-h-voffset/2],...
                'XAxisLocation', 'Top')
    
    ha = axes('parent', hf,...
              'Position',[rect_o(1) voffset rect_o(3) h],...
              chromaxes_props{:},  'Visible', 'off');
    orientation = 2;
else
    if isempty(hf_old)
        hf = figure('Color',[1 1 1],...
                    'Tag','Bioinfo:Karyotype');
        ha = axes('parent', hf,...
                  chromaxes_props{:});
        pos = get(ha, 'Position');
        
        set(ha,  'Position', [0 0 1 pos(2)+pos(4)]) 
        axis off;
    else
        %== Deactivate all figure UI modes
        activateuimode(hf_old, '')
        hf = hf_old;
    end
end

appdata = localGetAppData(hf);

appdata.orientation = orientation; % 1 is vertical and 2 is horizontal
appdata.chromnum = chromnum;
if unit == 1
    appdata.unit = 1;
    appdata.unitstr = '';
elseif unit == 2
    appdata.unit = 1000;
    appdata.unitstr = '(kb)';
elseif unit == 3
    appdata.unit = 10^6;
    appdata.unitstr = '(mb)';
end

if ~isempty(hf_old)
    delete(appdata.plotgrp)
    appdata.plotgrp = hggroup;
    if appdata.showcnv && ~isempty(appdata.cnvlines)
        delete(appdata.cnvlines)
		appdata.showcnv = false;
    end
else
    appdata.hf = hf;
    appdata.ha = ha;
    appdata.plotgrp = hggroup;
    appdata.cytobandstruct = cytobandstruct;
    appdata.addtoplot = addtoplot;
    % Edges
    if appdata.addtoplot
        appdata.vEDGE = 0;
        appdata.hEDGE = 0;
    else
        appdata.vEDGE = 0.05;
        appdata.hEDGE = 0.05;
    end
    appdata.showbandlabel = showbandlabels;

    % Chromosome ids
    appdata.giestains = label2GieStain(cytobandstruct.GieStains);
    appdata.chromosome = label2Chromosome(cytobandstruct.ChromLabels);
    appdata.chromStarts = double(cytobandstruct.BandStartBPs)/appdata.unit;
    appdata.chromEnds = double(cytobandstruct.BandEndBPs)/appdata.unit;
    appdata.bandLabels = cytobandstruct.BandLabels;
    [appdata.chromID, I] = unique(appdata.chromosome);
    
    % CNV
    if isempty(cnvstruct)
        appdata.showcnv = false;
    else
        appdata.showcnv = true;
        appdata.cnvstruct = cnvstruct;
        appdata.cnvchromosome = label2Chromosome(cnvstruct.Chromosome, cytobandstruct.ChromLabels);
        appdata.cnvType = cnvstruct.CNVType;
        appdata.cnvStarts = double(cnvstruct.Start)/appdata.unit;
        appdata.cnvEnds = double(cnvstruct.End)/appdata.unit;
    end

    appdata.bandlabels_maxWD =0;

    % Number of chromosomes - for example: 24 for Homo sapiens
    if chromnum == 0
        % SL CHANGE
        N = numel(appdata.chromID);
%         N = 23;
    else
        N = 1;
    end

    % Chromosome labels
    appdata.chromLabels = cytobandstruct.ChromLabels(I);

    if chromnum == 0
        [appdata.chromLText, appdata.chromLExt] = ...
            getChromLabelText(1:N, appdata.chromLabels);
        appdata.showbandlabel = false;
    else
        if ~appdata.addtoplot
            [appdata.chromLText, appdata.chromLExt] = ...
                getChromLabelText(chromnum, appdata.chromLabels);
        else
            appdata.chromLExt = 0;
        end
        appdata = getBandLabelText(appdata);
        appdata.bandlabel_markers = [];
    end

    % Determine the number of rows and columns
    if chromnum == 0
        if N/2 <= 15
            ncols = ceil(N/2);
        else
            % Draw 15 chromosomes in one row, get number of rows
            ncols = 15;
        end
        nrows = ceil(N/ncols);
    else
        nrows = 1;
        if appdata.addtoplot
            ncols = 1;
        else
            ncols = 9;
        end
    end

    % Create a matrix hold the start and end bp of the chromosomes
    if chromnum == 0
        appdata.chr_startend_bp = zeros(N,2);

        for i = 1:N
            idx = find(appdata.chromosome == i);
            appdata.chr_startend_bp(i, :) = [appdata.chromStarts(idx(1)),...
                appdata.chromEnds(idx(end))];
        end
    else
        idx = find(appdata.chromosome == chromnum);
        appdata.chr_startend_bp = [appdata.chromStarts(idx(1)),appdata.chromEnds(idx(end))];
    end

    appdata.max_bp = max(appdata.chr_startend_bp(:,2));
    min_bp = min(appdata.chr_startend_bp(:,1));
    
    if appdata.addtoplot
        set(ha_old, 'xlim', [min_bp appdata.max_bp])
    end

    appdata.chr_startend_bp = appdata.chr_startend_bp / appdata.max_bp;

    % The high and width of the area to draw one chromosome
    appdata.chrHeight = (1-2*appdata.vEDGE- appdata.chromLExt)/nrows;
    chr_WD = (1-2*appdata.hEDGE)/ncols;

    % Draw half width:
    if appdata.addtoplot
        appdata.draw_WD = (chr_WD - appdata.bandlabels_maxWD)/7;
        appdata.cnvGap = appdata.draw_WD/4;
    else
        appdata.draw_WD = chr_WD/6;
        appdata.cnvGap = appdata.draw_WD/2;
    end

    % Chromosome length in the plot
    appdata.chr_len = appdata.chr_startend_bp(:,2)' * appdata.chrHeight;

    %-----Plot Karyotype view-----------------------------------%
    % Draw the bottom lines of different row
    Y0 = (appdata.chrHeight*((nrows:-1:1)-1) + appdata.vEDGE)';
    X0 = (chr_WD*(0:ncols-1) + chr_WD/2 + appdata.hEDGE);

    X = repmat(X0, 2, nrows);
    if appdata.addtoplot
        xoffset =  0.16;
    else
        if appdata.chromnum == 0
            xoffset = 0;
        else
            xoffset = 0.35;
        end
    end
    appdata.X = X(:, 1:N)+ appdata.bandlabels_maxWD + xoffset;
    
    Y = repmat(Y0, 1, ncols);
    Y = reshape(Y', 1, ncols *nrows);
    Y = Y(1:N);

    appdata.Y = [Y+appdata.chr_len; Y] + appdata.chromLExt;
end

% Add chromosome labels
if ~appdata.addtoplot
	for i = 1:size(appdata.Y, 2)
		if appdata.orientation == 2
			set(appdata.chromLText(i),...
				'position', [1-appdata.Y(2,i)+appdata.chromLExt, appdata.X(1,i)],...
				'VerticalAlignment','Middle',...
				'visible', 'on');
		else
			set(appdata.chromLText(i), ...
				'position', [appdata.X(1,i), appdata.Y(2,i)-appdata.chromLExt/10],...
				'visible', 'on');
		end
	end
end

appdata = updateChromLim(appdata);
appdata = createContextMenus(appdata);
appdata = chromosomePatchs(appdata);

if appdata.showcnv
    appdata = drawCNVLines(appdata);
end

localSetAppData(hf,appdata);    

% Using datacursormode object
appdata.dcmobj = datacursormode(hf);
appdata.dcmUpdateFcnMainAxes = get(appdata.dcmobj,'UpdateFcn');
set(appdata.dcmobj, 'UpdateFcn', @updateBandDataCursor,...
                    'Enable','on')
hdcmMenu = get(appdata.dcmobj,'UIContextMenu'); % Get its context menu.
delete(findobj(hdcmMenu,'Tag','DataCursorEditText'))
delete(findobj(hdcmMenu,'Tag','DataCursorSelectText'));
set(appdata.dcmobj,'Enable','off') %starts by default off

% Using zoom mode object
if appdata.addtoplot
    appdata.zmobj = zoom(hf);
    setAllowAxesZoom(appdata.zmobj, appdata.ha, false)
    set(appdata.zmobj, 'ActionPostCallback', @updateAfterZoom)
    appdata.panobj = pan(hf);
    setAllowAxesPan(appdata.panobj,appdata.ha,false)
    setAxesPanMotion(appdata.panobj,ha_old,'Vertical')
    
	boxratio = 2.4;
	if appdata.showcnv
		boxratio = 2.8;
	end
		
    appdata.rect = rectangle('position', [1-appdata.Y(1),...
										  appdata.X(1)-1.2*appdata.draw_WD,...
										  appdata.chr_len,...
										  boxratio*appdata.draw_WD],...
                             'EdgeColor', [1 0.5 0.2],...
                             'Linewidth', 1.5,...
                             'Visible', 'off');
end

% Save state
appdata.datatip = getDataTipText();
localSetAppData(hf,appdata);
displayBandLabels([], hf);
set(hf,'WindowButtonMotionFcn',@localWindowButtonMotion);

if nargout > 0
    hfchrom = hf;
end
end % chromsomeplot function

%--------------------Helper functions-----------------------%
function appdata = chromosomePatchs(appdata)
% Create stain patches

x_ini = appdata.X(1,:);
y_ini = appdata.Y(1,:);

N = size(x_ini,2);
chr_acent_bp = zeros(3,N);
wd = appdata.draw_WD; % chromosome width
ht =  appdata.chr_len; % chromosome height

for i = 1:N
    x0 = x_ini(i);
    y0 = y_ini(i);

    blx = [];
    bly = [];
    if appdata.chromnum == 0
        idx = find(appdata.chromosome == i);
        chr_label = appdata.chromLabels{i};
    else
        idx = find(appdata.chromosome == appdata.chromnum);

        blx = zeros(2, numel(idx));
        bly = zeros(2, numel(idx));
        chr_label = appdata.chromLabels{appdata.chromnum};
    end 
    chr_giestains =appdata.giestains(idx);
    chr_start = appdata.chromStarts(idx);
    chr_end = appdata.chromEnds(idx);
    chr_band = appdata.bandLabels(idx);
    
    % Find centromere
    p_centr = [];
    cidx = find(chr_giestains == 8);
    if ~isempty(cidx)        
        chr_acent_bp(1:2, i) = chr_start(cidx);
        chr_acent_bp(3, i) = chr_end(cidx(end));
        
        p_centr = appdata.chrHeight*chr_acent_bp(:,i)/appdata.max_bp;
        centr_patchs = bandShape(x0, y0, wd, ht(i), p_centr, 8, appdata.orientation);
        for c = 1:numel(cidx)
            band.chrom = chr_label;
            band.start = chr_start(cidx(c));
            band.end = chr_end(cidx(c));
            band.label = chr_band{cidx(c)};
            set(centr_patchs(c), 'Userdata', band, 'parent', appdata.plotgrp);
        end
        
        if appdata.chromnum ~= 0
            for c=1:numel(cidx)
                [blx(:, cidx(c)), bly(:, cidx(c))] = setBandLabels(...
                    centr_patchs(c),appdata.bandlabels_text(cidx(c)),...
                    appdata.cnvGap, appdata.orientation);
            end
        end
        set(centr_patchs, 'UIContextMenu', appdata.contextmenu)
    end
    
    % SL CHANGE
    % Find stain bands
%     n = length(find(chr_giestains >0 & chr_giestains < 8));
%     gpatch = zeros(n, 1);
%     ncount = 1;
%     for g = 1:7
%        cidx = find(chr_giestains == g);
%        if ~isempty(cidx)
%            yc =appdata.chrHeight* [chr_start(cidx)'; chr_end(cidx)']/appdata.max_bp;
%            for c = 1:length(cidx)
%                gpatch(ncount) = bandShape(x0, y0, wd, ht(i), yc(:, c), g, appdata.orientation);
%                band.chrom = chr_label;
%                band.start = chr_start(cidx(c));
%                band.end = chr_end(cidx(c));
%                band.label = chr_band{cidx(c)};
%                set(gpatch(ncount), 'Userdata', band);
%                if appdata.chromnum ~= 0
%                    [blx(:, cidx(c)), bly(:, cidx(c))] = setBandLabels(...
%                     gpatch(ncount),appdata.bandlabels_text(cidx(c)),...
%                     appdata.cnvGap, appdata.orientation); %#ok<AGROW>
%                end
%                ncount = ncount+1;
%            end
%        end
%     end
%     set(gpatch, 'UIContextMenu', appdata.contextmenu,  'parent', appdata.plotgrp)

    % Find stalk
    cidx = find(chr_giestains == 9);
    p_stalks = [];
    if ~isempty(cidx)
        p_stalks =appdata.chrHeight* [chr_start(cidx)'; chr_end(cidx)']/appdata.max_bp;
        for c = 1:length(cidx)
            stalkpatch = bandShape(x0, y0, wd, ht(i), p_stalks(:, c), 9, appdata.orientation);
            band.chrom = chr_label;
            band.start = chr_start(cidx(c));
            band.end = chr_end(cidx(c));
            band.label = chr_band{cidx(c)};
            set(stalkpatch, 'Userdata', band);
            if appdata.chromnum ~= 0
                [blx(:, cidx(c)), bly(:, cidx(c))] = setBandLabels(...
                    stalkpatch,appdata.bandlabels_text(cidx(c)),...
                    appdata.cnvGap, appdata.orientation); %#ok<AGROW>
            end
            set(stalkpatch, 'UIContextMenu', appdata.contextmenu,  'parent', appdata.plotgrp);
        end        
    end

    appdata.bandlabel_markers = line(blx, bly, 'color', 'k', 'visible', 'off');
    set(appdata.bandlabel_markers, 'parent', appdata.plotgrp)
    % Draw outline
	if ht(i) >= 2*wd
		chromOutline(x0, y0, wd, ht(i), p_stalks, p_centr, appdata );
	end
end
end % chromosomePatchs function
%--------------------------------------------------------------------%
function chromOutline(xo, yo, wd, ht, p_stalk, p_centr, appdata)
% Return the outline line object of a chromosome
% p_stalk - stalk positions
% p_centr - centromere positions

cdir = appdata.orientation;
lineprops = {'color', [0 0 0],...
               'linewidth', 0.5,...
               'parent', appdata.plotgrp};
 
% Draw the chromosome ends
[xc, yc] = curveCoords(wd, [1 0]');

xc = [0; xc];
x = repmat(xc, 1, 4);
x(:,[2, 4]) = flipud(x(:,[2, 4]));
x(:,1:2) = xo+x(:,1:2);
x(:,3:4) = xo-x(:,3:4);

yc = [0; yc];
y = repmat(yc,1,4);
y(:,[2,4]) = flipud(y(:,[2,4]));
y(:, [1,4]) = yo-y(:, [1,4]);
y(:, [2,3]) = yo - ht +y(:, [2,3]);

if cdir ==2
    [xe, ye] = flipXY(x,y);
    line(xe, ye, lineprops{:});
else
    line(x, y, lineprops{:});
end

% draw the side lines
if isempty(p_stalk) && isempty(p_centr)
    xs = [x(end, 1)  x(end, 3);...
        x(1, 2) x(1,4)];
    ys = [y(end, 1)  y(end, 3);...
        y(1, 2) y(1,4)];

    if cdir ==2
        [xs, ys] = flipXY(xs,ys);
    end
    line(xs, ys, lineprops{:});
    return;
end

if isempty(p_stalk)
    p_centr = yo - p_centr;
    xs = [x(end, 1), x(1, 2), x(end, 3), x(1,4);...
        x(end, 1), x(1, 2), x(end, 3), x(1,4)];
    ys = [y(end, 1), p_centr(3), y(end, 3), p_centr(1);...
        p_centr(1), y(1, 2), p_centr(3), y(1,4)];

    if cdir ==2
        [xs, ys] = flipXY(xs,ys);
    end
    line(xs, ys, lineprops{:});
    return;
end

p = p_stalk;
if ~isempty(p_centr) % add it to p_stalk
    idx = find(p_stalk(2,:) < p_centr(1), 1, 'last');
    p = [p_stalk(:, 1:idx), p_centr([1,3]), p_stalk(:, idx+1:end)];
end

p = yo - p;
yss = [[1 y(end, 1)]', p, [y(1,2) 1]'];

% Flip yss and shift to right on first row
yss = flipud(yss);
yss_1 = yss(1,:);
yss(1,:) = yss_1([end 1:end-1]);
yss = yss(:, 2:end);

xss = [ones(size(yss)) * x(end,1), ones(size(yss)) * x(1,4)];
yss = [yss, fliplr(yss)];

if cdir ==2
    [xss, yss] = flipXY(xss,yss);
end
line(xss, yss, lineprops{:});
return;
end % chromOutline function
%-------------------------------------------------------------------------%
function shapeobj = bandShape(xo, yo, wd, ht, p, gtype, cdir)
% This function return a patch object of certain shape for a band. The
% shapes have no edge outline. 
% xo, yo - the origin of the chromosome
% p - two or three element column vector of y start and end points
% wd - the width of the chromosome
% ht - height of the chromosome
% gtype - giestain type: 
%               gneg     - 1
%               gpos25   - 2
%               gpos50   - 3
%               gpos75   - 4
%               gpos100  - 5
%               gpos     - 6
%               gvar     - 7
%               acen     - 8
%               stalk    - 9


% Get the centromere red triangles
if gtype == 8
    coord1 = @(x,w) repmat([x+wd, x, x-wd]', 1, 2);
    coord2 = @(y,pt) y - [pt(1) pt(3); pt(2) pt(2); pt(1) pt(3)];
    x = coord1(xo, wd);
    y = coord2(yo, p);
    
    if cdir == 2
        [x, y] = flipXY(x,y);
    end

    c = [0 0 0]; % red % SL CHANGE TO BLACK

    shapeobj = [patch(x(:,1), y(:,1), c, 'edgecolor', 'none'),...
        patch(x(:,2), y(:,2), c, 'edgecolor', 'none')];
    return;
end

if gtype == 9 % stalk
    delta = (p(1) - p(2))/5;    
    coord1 = @(x, w)[x+w; x-w];
    coord2 = @(y, h)[y - h; yo-flipud(h)];
    
    h = [p(1) p(1)-delta p(2)+delta p(2)]';
    w = [wd wd/2 wd/2 wd]';
    x = coord1(xo, w);
    y = coord2(yo, h);
    
    if cdir == 2
       [x, y] = flipXY(x,y);
    end

    c =[0.9 0.9 0.9];
    shapeobj = patch(x, y, c, 'edgecolor', 'none');
    return;
end

if any(gtype==1:7)
    % Determine the curve bounds
    bounds = @(y, w, h) y -h + w/2;
    p_bound = bounds(yo, -wd, 0);
    q_bound =  bounds(yo, wd, ht);
    
    isonparm = @(x, pt, bd) x - pt(1) > bd;
    isonqarm = @(x, pt, bd) x - pt(2) < bd;
    
    coord1 = @(x, c) [x + c; x - flipud(c)];
    coord2 = @(y, c) [y + c; y + flipud(c)];
    
    % Get stain bands
        if isonparm(yo, p, p_bound)
            [cv1,cv2] = curveCoords(wd, wd/2 - p);
            if isonqarm(yo, p, p_bound)
                cv1 = [cv1; wd];
                cv2 = [cv2; p(2)];
            end

            x = coord1(xo, cv1);
            y = coord2(yo, -cv2);
        elseif isonqarm(yo, p, q_bound)
            [cv1,cv2] = curveCoords(wd, q_bound - yo+[p(2) p(1)]');
            cv1 = flipud(cv1);
            cv2 = flipud(cv2);

            if isonparm(yo, p, q_bound)
                cv1 = [wd; cv1];
                cv2 = [ht-p(1); cv2];
            end
            x = coord1(xo, cv1);
            y = coord2(yo-ht, cv2);
        else
            x = coord1(xo, ones(2,1)*wd);
            y = coord2(yo, -p);
        end

    if cdir == 2
        [x, y] = flipXY(x,y);
    end
    
    switch gtype
        case 1 %- gneg
            c = [1 1 1];
        case 2 %- gpos25
            c = [0.8 0.8 0.8];
        case 3 %- gpos50
            c = [0.5 0.5 0.5];
        case 4 %- gpos75
            c = [0.3 0.3 0.3];
        case {5,6} %- gpos100 and gpos
            c = [0 0 0];
        case 7 % gvar
            c = [0.8 0.9  1];
    end
    
    shapeobj = patch(x, y, c, 'edgecolor', 'none');
    return;  
end
end % bandShape function
%-------------------------------------------%
function [cx, cy] = curveCoords(wd, delta)
% Return the x,y coordinates of an arc of 90 degrees.
% wd - width of chromosome
% corner - 1,2,3 and 4 for the four coners

np = 5;
r = wd/2;
angle1 = asin(min(delta(1)/r, 1));
angle2 = asin(max(delta(2)/r, 0));

theta = linspace(angle1,angle2,np);
rho = ones(size(theta))*r;
[cx,cy] = pol2cart(theta, rho);
cx = (cx + wd-r)';
cy = (r - cy)'; 
end % curveCoords function
%---------------------------------------------%
function appdata = getBandLabelText(appdata)
% Get band labels text information of a chromosome

idx = find(appdata.chromosome == appdata.chromnum);
N = numel(idx);
appdata.bandlabels_text = zeros(N, 1);
appdata.bandlabels_extents = zeros(N, 1);

fontsize = 7;
if strncmpi(computer, 'MAC', 3)
    fontsize = 8;
end

bandLabels = appdata.bandLabels(idx);
for i = 1:numel(idx)
    appdata.bandlabels_text(i) = text(0,1,1,bandLabels{i},...
        'Color', [0 0 0],...
        'Margin', 1,...
        'VerticalAlignment','Middle',...
        'HorizontalAlignment','right',...
        'Clipping','on',...
        'Interpreter','none',...
        'Fontsize', fontsize,...
        'Visible','off');
    appdata.bandlabels_extents(i)= get(appdata.bandlabels_text(i), 'Extent')*[0;0;1;0];
end
appdata.bandlabels_maxWD = max(appdata.bandlabels_extents);
end % getBandLabelText function
%-------------------------------------------------
function [labels, maxext] = getChromLabelText(chromNums, chromLabels)
% chromNums - The number of chromosomes to create label for
% chromLabels - labals for chromsome

N = numel(chromNums);
labels = zeros(N,1);
label_extents = zeros(N,1);
for i = 1:N
    labels(i) = text(0,1,1,chromLabels{chromNums(i)},...
                'Tag','Chromlabel',...
                'Color', [0 0 0],...
                'Margin', 1,...
                'VerticalAlignment','top',...
                'HorizontalAlignment','Center',...
                'Clipping','on',...
                'Interpreter','none',...
                'Visible','off');
    label_extents(i) = get(labels(i), 'Extent')*[0;0;0;1];
end
maxext = max(label_extents);
end  % getChromLabelText function
%-------------------------------------------------
function datatip = getDataTipText()
fs = 8;
if strncmpi(computer, 'MAC', 3)
    fs = 9;
end
datatip = text(0,1,1,'k',...
        'Tag','chromdatatip',...
        'BackgroundColor',[1 1 .93],...
        'Color', [0 0 0],...
        'EdgeColor', [0.8 0.8 0.8],...
        'VerticalAlignment', 'Top',...
        'Clipping','off',...
        'Visible','off',...
        'Fontsize',fs,...
        'Interpreter','none');
end % getDataTipText function
%---------------------------------------------%
function appdata = drawCNVLines(appdata)
% Draw CNV lines

xo = appdata.X(1,:);
yo = appdata.Y(1,:);

if appdata.chromnum == 0
    idx = true(numel(appdata.cnvchromosome),1);
else
    idx = appdata.cnvchromosome == appdata.chromnum;
end

cnvchromosome = appdata.cnvchromosome(idx);
cnvstart = appdata.chrHeight*appdata.cnvStarts(idx)/appdata.max_bp;
cnvend = appdata.chrHeight*appdata.cnvEnds(idx)/appdata.max_bp;
cnvtype = appdata.cnvType(idx);
[ucntype, typeidx] = unique(cnvtype);
appdata.cnvlines = zeros(numel(cnvchromosome),1);
cnvlegh = zeros(1, numel(ucntype));
cnvTypes = {'Loss', 'Gain'};

for i=1:numel(cnvchromosome)
   chromnum = cnvchromosome(i);
   
   if appdata.chromnum == 0
       x0 = xo(chromnum);
       y0 = yo(chromnum);
   else
       x0 = xo;
       y0 = yo;
   end
   marker.start = appdata.cnvStarts(idx(i));
   marker.end = appdata.cnvEnds(idx(i));
   
   y = [y0 - cnvstart(i); y0 - cnvend(i)];
   if cnvtype(i) == 1 % Loss
       x = (x0 - appdata.cnvGap - appdata.draw_WD) * [1;1];
       c = [ 1 0 0]; % red
       marker.type = 'Loss';
   elseif cnvtype(i) == 2 % Gain
       x = (x0 + appdata.cnvGap + appdata.draw_WD) * [1;1];
       c = [0.2 0.8 .2];
       marker.type = 'Gain';
   end

   if appdata.orientation == 2
       [x, y] = flipXY(x,y);
   end
   appdata.cnvlines(i) = line(x, y, 'Color', c,...
                                    'Linewidth', 2,...
                                    'Tag', 'cnvs',...
                                    'Userdata', marker,...
                                    'Visible', 'on');
   hasbehavior(appdata.cnvlines(i), 'legend', false);
end

for i = 1:length(cnvlegh)
    cnvlegh(i) = appdata.cnvlines(typeidx(i));
    hasbehavior(cnvlegh(i), 'legend', true);
    set(cnvlegh(i),'DisplayName',cnvTypes{ucntype(i)});
end
end % drawCNVLines function
%---------------------------------------------%
function localWindowButtonMotion(hfig, varargin)
% Callback function activated when moving over the axes, checks location of
% the mouse and puts datatip if over an active node.
appdata = localGetAppData(hfig);

dcmobj = datacursormode(hfig); %get handle to the DataCursor in current figure

ison = get(dcmobj, 'Enable');
if strcmpi(ison, 'on')
    return;
end
% set a virtual grid to get the point
xThres=appdata.chr_xlim;
yThres=appdata.chr_ylim;
cp = get(appdata.ha,'CurrentPoint');
xPos = cp(1,1); yPos = cp(1,2);

hp = xPos>xThres(1,:) & xPos<xThres(2,:) & ...
     yPos>yThres(1,:) & yPos<yThres(2,:);
hp = find (hp);

online = [];
if appdata.showcnv && isempty(hp)
    if numel(appdata.cnvlines) > 1
        xvals = (cell2mat(get(appdata.cnvlines, 'XData')))';
        yvals = (cell2mat(get(appdata.cnvlines, 'YData')))';
    else
        xvals = (get(appdata.cnvlines, 'XData'))';
        yvals = (get(appdata.cnvlines, 'YData'))';
    end
    
    if isempty(xvals) || isempty(yvals)
        online = [];
    else
        xrange = get(appdata.ha,'Xlim');
        yrange = get(appdata.ha,'Ylim');

        fuzzx = 0.01 * (xrange(2) - xrange(1));
        fuzzy = 0.01 * (yrange(2) - yrange(1));
        online = find(((yPos > yvals(2,:) - fuzzy & yPos < yvals(1,:) + fuzzy) &...
            (xPos > xvals(1,:) - fuzzx & xPos < xvals(2,:) + fuzzx)), 1);
    end
end

if isempty(hp) && isempty(online)
    set(appdata.datatip, 'visible', 'off')
else
    if appdata.unit == 1
        unitfmt = '%d';
    else
        unitfmt = '%0.1f';
    end

    if ~isempty(hp)
        if appdata.chromnum == 0
            chromnum = hp;
        else
            chromnum = appdata.chromnum;
        end
        
        idx = find(appdata.chromosome == chromnum);
        labels = {sprintf('Chromosome %d', chromnum);...
            'pter - qter';...
            sprintf(['Start %s: ' unitfmt], appdata.unitstr, appdata.chromStarts(idx(1)));...
            sprintf(['End %s: ' unitfmt], appdata.unitstr, appdata.chromEnds(idx(end)))};
    else
        if appdata.cnvType(online) == 1
            cnvtype = 'Loss';
        else
            cnvtype = 'Gain';
        end
   
        if appdata.unit == 1
            bpstr =[num2str(appdata.cnvStarts(online)),...
                    '-', num2str(appdata.cnvEnds(online))];
        else
            bpstr =[num2str(appdata.cnvStarts(online), '%.1f'),...
                    '-' num2str(appdata.cnvEnds(online), '%.1f'), appdata.unitstr];
        end
        
        labels = {sprintf('CNV on Chromosome %d', appdata.cnvchromosome(online));...
                           cnvtype; bpstr};
        if appdata.chromnum == 0
            hp = appdata.cnvchromosome(online);
        else
            hp = 1;
        end
    end
    set(appdata.datatip, 'String', char(labels));
    extents = get(appdata.datatip, 'Extent');
    xlim = get(appdata.ha, 'xlim');

    if appdata.orientation == 1
        y = yPos;
        if yPos - yThres(1, hp)< extents(4)
            y = yThres(1, hp)+extents(4);
        end

        x = xThres(2, hp)+appdata.draw_WD;
        if x + extents(3) > xlim(2)
            x = xThres(1, hp) - appdata.draw_WD/2 - extents(3);
        end
        vertalign = 'Top';
    else
        x = xPos;
        if xPos + extents(3) > xThres(2, hp)
            x = xThres(2, hp)- extents(3);
        end

        y = yThres(2, hp)+appdata.draw_WD;
        vertalign = 'Bottom';
    end

     set(appdata.datatip,...
        'position', [x, y],...
        'VerticalAlignment', vertalign,...
        'Visible', 'on');
end
end % localWindowButtonMotion( function
%---------------------------------------------------%
function appdata = createContextMenus(appdata)
% helper function to set context menus for the chromosome groups
% 
% hcm   contains handle to the chromosome context menu to
%        used to reactivate it

% context menu for chromosomes
hcm = uicontextmenu('Callback', @enableContextMenu);

if appdata.chromnum == 0
    item = uimenu(hcm,'Label','Display in New Figure',...
                      'Tag', 'newfigcontrol');
    uimenu(item,'Label','Vertical','Callback',{@displaySingleChrom,1});
    uimenu(item,'Label','Horizontal','Callback',{@displaySingleChrom,2});
else
    uimenu(hcm,'Label','Show G-band Labels',...
       'Checked', 'off',...
       'Tag', 'bandlabelcontrol',...
       'Callback',@displayBandLabels);
   if ~appdata.addtoplot
       item = uimenu(hcm,'Label','Display',...
           'Tag', 'scdisplaycontrol');
       checkstate = {'on', 'off'};
       idx = appdata.orientation == [1 2];

       uimenu(item,'Label','Vertical',...
           'Checked', checkstate{idx},...
           'Callback', {@updateOrientation,1});
       uimenu(item,'Label','Horizontal',...
           'Checked', checkstate{~idx},...
           'Callback',{@updateOrientation,2});
   end
end

appdata.contextmenu = hcm;
end % createContextMenus function
%-----------------------------------------------------------%
function enableContextMenu(h, varargin)
appdata = localGetAppData(gcbf);

set(appdata.datatip, 'visible', 'off')
hc = get(h,'Children');

if strcmpi(hc, 'bandlabelcontrol') && appdata.showbandlabel
    set(hc, 'Label', 'Hide G-band Labels');
    set(appdata.bandlabels_text,'visible', 'on')
    set(appdata.bandlabel_markers, 'visible', 'on')
end
set(hc,'Enable','on')
end % enableContextMenu function
%-----------------------------------------------------------%
function displayBandLabels(h, varargin)
% Show bandlabels in single chromosome mode

if isempty(h)
    hfig = varargin{1};
else
    hfig = gcbf;
end

appdata = localGetAppData(hfig);

if appdata.chromnum == 0
    return;
end

if ~isempty(h)
    isshown = strncmp(get(h, 'Label'), 'Hide', 4);
    if isshown
        appdata.showbandlabel = false;
        set(h, 'Label', 'Show G-band Labels');
    else
        appdata.showbandlabel = true;
        set(h, 'Label', 'Hide G-band Labels');
    end
    localSetAppData(hfig, appdata);
else
    h = findobj(get(appdata.contextmenu,'Children'),'Tag','bandlabelcontrol');
    if appdata.showbandlabel
        set(h, 'Label', 'Hide G-band Labels');
    else
        set(h, 'Label', 'Show G-band Labels');
    end
end

if appdata.showbandlabel
    set(appdata.bandlabels_text,'visible', 'on')
    set(appdata.bandlabel_markers, 'visible', 'on')
else
    set(appdata.bandlabels_text,'visible', 'off')
    set(appdata.bandlabel_markers, 'visible', 'off')
end

end % displayBandLabels function
%---------------------------------------------------------%
function displaySingleChrom(h, ent, varargin)  %#ok<INUSL>
% Display a single chromosome in a new figure
hf_old = gcbf;
appdata = localGetAppData(hf_old);

cdir = varargin{1};

% Retrieve the chromosome number clicked on
datatip = get(appdata.datatip, 'String');
chromnum = strtrim(datatip(1,:));
chromnum = str2double(chromnum(12:end));

if appdata.addtoplot
    chromosomeplot(appdata.cytobandstruct, chromnum,...
        'orientation', cdir);
else
	if appdata.showcnv
		chromosomeplot(appdata.cytobandstruct, chromnum,...
			'cnv', appdata.cnvstruct,...
			'orientation', cdir);
	else
		chromosomeplot(appdata.cytobandstruct, chromnum,...
			'orientation', cdir);
	end
end
end % displaySingleChrom function
%---------------------------------------------------%
function appdata = updateChromLim(appdata)
% Update chromsome xlim and ylim according to orientation

if appdata.orientation == 1
    appdata.chr_xlim = [appdata.X(1,:)-appdata.draw_WD; appdata.X(1,:)+appdata.draw_WD];
    appdata.chr_ylim = [appdata.Y(end,:); appdata.Y(1,:)];
else
    appdata.chr_xlim = [appdata.Y(end,:); appdata.Y(1,:)];
    appdata.chr_ylim = [appdata.X(1,:)-appdata.draw_WD; appdata.X(1,:)+ appdata.draw_WD];
end
end % updateChromLim function
%---------------------------------------------%
function updateOrientation(h, ent, varargin) %#ok<INUSL>
% Display a single chromosome in a new figure
% Display a single chromosome in a new figure
hf_old = gcbf;
appdata = localGetAppData(hf_old);

cdir = varargin{1};

% Retrieve the chromosome number clicked on
if appdata.showcnv
	chromosomeplot(appdata.cytobandstruct, appdata.chromnum,...
		'cnv', appdata.cnvstruct,...
		'orientation', cdir,...
		'figure', hf_old);
else
	chromosomeplot(appdata.cytobandstruct, appdata.chromnum,...
		'orientation', cdir,...
		'figure', hf_old);
end
end % updateOrientation function
% -----------------------------------------%
function localSetAppData(hfig,appdata)
setappdata(hfig,'IdeogramPlot',appdata);
end % localSetAppData function
%----------------------------------------------%
function [appdata] = localGetAppData(hfig)
if isappdata(hfig,'IdeogramPlot')
    appdata = getappdata(hfig,'IdeogramPlot');
else
    appdata = guihandles(hfig);
    appdata.organism = [];
    appdata.numchrom = [];
end
end % localGetAppData function
%-------------------------------------------------%
function chrs = label2Chromosome(chromLabels, cytoChromLabels)
% Convert chromosome label '2' to a number 2
if nargin < 2
    if isnumeric(chromLabels)
        chrs = chromLabels(:);
    else
        chrs = str2double(chromLabels);
        chrs(strcmpi('X', chromLabels)) = max(chrs)+1;
        chrs(strcmpi('Y', chromLabels)) = max(chrs)+1;
    end
else
    cytoChromLabels = unique(cytoChromLabels);
    chrs = zeros(length(chromLabels), 1);
    for i = 1:length(chromLabels)        
        % Convert to number
        if isnumeric(chromLabels) 
            chrLabel = chromLabels(i);
        else
            chrLabel = strtrim(chromLabels(i, :));
            tmp = str2double(chromLabels(i, :)); %ok
            if ~isnan(tmp)
                chrLabel = tmp;
            end
        end
        
        chrs(i) = validateChromNumber(cytoChromLabels, chrLabel);
        if chrs(i) < 0
            if strcmpi(get(gcf, 'Name'), '')
                close(gcf)
            end
            error(message('bioinfo:chromosomeplot:CNVChromNumNotValid', num2str( chromLabels( i ) )));
        end
    end
end
chrs = int8(chrs);
end % label2Chromosome function
% ------------------------------------------------
function chrnum = validateChromNumber(chromLabels, chrnum)
% Validate chromnum input variable
if isempty(chrnum)
    chrnum = 0;
else
    N = numel(chromLabels);
    if ischar(chrnum)
        if isempty(find(strcmpi(chrnum, chromLabels), 1))
            chrnum = -1;
        else
            chrnum = (N-2) + find(strcmpi(chrnum, {'X', 'Y'}));
        end
    elseif chrnum > N
        chrnum = -1;
    end
end
end % validateChromNumber function
% ------------------------------------------------%
function giestain = label2GieStain(gielabel)
% Convert giestain  label 'gneg' to a number 1
% gtype - giestain type: 
%               gneg     - 1
%               gpos25   - 2
%               gpos50   - 3
%               gpos75   - 4
%               gpos100  - 5
%               gpos     - 6
%               gvar     - 7
%               acen     - 8
%               stalk    - 9
N = numel(gielabel);
giestain = zeros(size(gielabel));
for i=1:N
	x = gielabel{i};
    if strcmpi(x, 'gneg')
        giestain(i) = 1;
    elseif strncmp(x, 'gpos', 4)
        n = str2double(strtok(x, 'gpos'));
        if n < 50  %strcmpi(x, 'gpos25')
            giestain(i) = 2;
        elseif n < 75 && n >= 50% strcmpi(x, 'gpos50')
            giestain(i) = 3;
        elseif n < 100 && n >= 50% strcmpi(x, 'gpos75')
            giestain(i) = 4;
        elseif n>= 100 % strcmpi(x, 'gpos100')
            giestain(i) = 5;
        elseif isnan(n) % no gpos number
            giestain(i) = 6;
        end
    elseif strcmpi(x, 'gvar')
        giestain(i) = 7;
    elseif strcmpi(x, 'acen')
        giestain(i) = 8;
    elseif strcmpi(x, 'stalk')
        giestain(i) = 9;
    end
end
giestain = int8(giestain);
end % label2GieStain function
%------------------------------------------------%
function [blx, bly, mx, my]= getBandLabelXY(bandpatch, cnvgap, cdir)
% Return bandlabel blx,bly coordinates and marker mx, my
% wd - chromosome half width

py = get(bandpatch, 'YData');
px = get(bandpatch, 'XData');

halfwidth = @(x)(max(x) - min(x))/2;
if cdir == 1
    wd = halfwidth(px);
else
    wd = halfwidth(py);
end
markerwd = wd*3/4;

if cdir == 1
    bly = getMidPoint(py);
    blx =  min(px) - markerwd - 1.5*cnvgap;
    mx = [blx+markerwd; blx];
    my = [bly;bly];
else
    blx = getMidPoint(px);
    bly =  min(py) - markerwd - 1.5*cnvgap;
    mx = [blx; blx];
    my = [bly; bly+markerwd];
end
end % getBandLabelXY function
%--------------------------------------------------------%
function mpt = getMidPoint(pt)
% Return the midpoint of a patch
if size(pt, 1) ==3 % centromere
    mpt = pt(1) + [1 0]*diff(pt)/2;
else
    nl = length(pt)/2;
    mpt = (pt(1, :) + pt(nl,1))/2;
end
end % getMidPoint function
%----------------------------------------%
function [mx, my] = setBandLabels(pobj, tobj, cnvgap, cdir)
% Set band label text object properties, and return labl marker position
% pobj - patch object
% tobj - label text obj
% x,y - positions
% cidr - orientations

setlabel = @(t, x, y, theta) set(t,...
                    'position', [x, y],...
                    'rotation', theta);
                
[lx, ly, mx, my]= getBandLabelXY(pobj,cnvgap,cdir); 

if cdir == 1
    setlabel(tobj, lx,ly,0);
else
    setlabel(tobj, lx,ly,90);
end

end % setBandLabels function
%-----------------------------------------------------%
function updateAfterZoom(obj, evt)  %#ok<INUSD>
% After zoom into the plot. Show update the rectangle box to indicate the
% region

ax = evt.Axes;
appdata = localGetAppData(obj);
xlim = get(ax, 'Xlim')/appdata.max_bp ;
pos = get(appdata.rect, 'Position');
set(appdata.rect, 'Position', [xlim(1) pos(2) xlim(2)-xlim(1) pos(4)])

if(diff(xlim) == 1)
    set(appdata.rect, 'Visible', 'off')
else
    set(appdata.rect, 'Visible', 'on')
end
end % updateAfterZoom function
%----------------------------------------%
function h = findParentAxes(h)
while ~strcmp(get(h,'Type'),'axes') && h~=0
    h = get(h,'Parent');
end
end
%----------------------------------%
function txt = updateBandDataCursor(empt, event_obj) %#ok<INUSL>
% Update data cursor display
target = get(event_obj, 'Target');
plotaxes = findParentAxes(target);

% Check if data point is not on the axes for chromosome plot or cnvs
% annotations
if ~strcmpi(get(plotaxes, 'Tag'), 'chromplotaxes')
    % honors previous datacursor function
    appdata = localGetAppData(get(findobj(gcbf,'Tag','chromplotaxes'),'parent'));
    if isempty(appdata.dcmUpdateFcnMainAxes)
         pos = get(event_obj, 'position');
         txt = {['X=',num2str(pos(1))], ['Y=',num2str(pos(2))]};
    else
         txt = feval(appdata.dcmUpdateFcnMainAxes,empt,event_obj);
    end
    return;
end

if strcmpi(get(target,'Type'), 'line') && strcmpi(get(target,'Tag'), 'cnvs')
    iscnvflag = true;
    %plotaxes = get(target, 'parent');
else
    iscnvflag = false;
    %plotaxes = get(get(target, 'parent'), 'parent');
end
  
if ~strcmpi(get(target, 'Type'), 'patch') && ~iscnvflag
    txt = {''};
    return;
end

hfig = get(plotaxes, 'parent');
appdata = localGetAppData(hfig);

if appdata.unit ==1
    unitfmt = '%d';
else
    unitfmt = '%.1f';
end

if iscnvflag
    marker = get(target, 'userdata');
    if isempty(marker)
        txt = {''};
    else
        txt = {marker.type,...
            [num2str(marker.start, unitfmt), '-',...
			num2str(marker.end, unitfmt), appdata.unitstr]};
    end
    
else
    band = get(target, 'userdata');
    if isempty(band)
        txt = {''};
        return;
    end

    txt = {['Chromosome ' band.chrom],...
        band.label,...
        [num2str(band.start, unitfmt), '-', num2str(band.end, unitfmt), appdata.unitstr]};
end
end % updateBandDataCursor function    
%-------------------------------------------------
function [x, y] = flipXY(x,y)
tmpx = x;
tmpy = y;

x = 1-tmpy;
y = tmpx;
end % flipXY function
%-----------------------------------------------------------%
function valid = checkCytoBandStruct(cbStruct)
% Checks for cytoband structure field names

valid =  isfield(cbStruct,'ChromLabels')&& ...
         isfield(cbStruct,'BandStartBPs') && ...
         isfield(cbStruct,'BandEndBPs') && ...
         isfield(cbStruct,'BandLabels') && ...
         isfield(cbStruct,'GieStains');
end % checkCytoBandStruct function
%-----------------------------------------------------------%
function valid = checkCNVStruct(cbStruct)
% Checks for CNV structure field names

valid =  isfield(cbStruct,'Chromosome')&& ...
         isfield(cbStruct,'CNVType') && ...
         isfield(cbStruct,'Start') && ...
         isfield(cbStruct,'End');
     
if valid
    if isnumeric(cbStruct.Chromosome)
        nchrs = numel(cbStruct.Chromosome);
    elseif ischar(cbStruct.Chromosome)
        nchrs = size(cbStruct.Chromosome, 1);
    else
        error(message('bioinfo:chromosomeplot:InvalidChromInput'));
    end
    ncnvs = [numel(cbStruct.CNVType), numel(cbStruct.Start), numel(cbStruct.End)];
    if any(ncnvs ~= nchrs)
        error(message('bioinfo:chromosomeplot:MismatchNumberCNVType', nchrs));
    end
end   
end % checkCNVStruct function
