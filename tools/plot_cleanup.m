function plot_cleanup(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modifies plots to make them look slightly better
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for vc = 1:2:length(varargin)
  switch(varargin{vc})
    case('AxisHandle')
      axishandle = varargin{vc+1};
    case('FigHandle')
      fighandle = varargin{vc+1};
    case('LineWidth')
      linewidth = varargin{vc+1};
		case('FontName')
			fontname = varargin{vc+1};
    case('FontSize')
      fontsize = varargin{vc+1};
		case('PlotBoxAspectRatio')
			pbaratio = varargin{vc+1};
		case('NoEqual')
			noequal = varargin{vc+1};
		case('Interpreter')
			interpreter = varargin{vc+1};
    otherwise
      warning('Unknown argument: %s',varargin{vc});
  end
end

if(~exist('axishandle','var')); axishandle = gca; end
if(~exist('fighandle','var')); fighandle = gcf; end
if(~exist('linewidth','var')); linewidth = 3; end
if(~exist('fontname','var')); fontname = 'Arial'; end
if(~exist('fontsize','var')); fontsize = 24; end
if(~exist('noequal','var')); noequal = false; end
if(~exist('pbaratio','var')); pbaratio = [1,1,1]; end
if(~exist('interpreter','var')); interpreter = 'tex'; end

set(fighandle,'defaulttextinterpreter',interpreter)

if(~noequal)
	set(axishandle,'Box','On','LineWidth',2,'PlotBoxAspectRatioMode',...
  'Manual','PlotBoxAspectRatio',pbaratio,'FontName',fontname,...
  'FontSize',fontsize,'TickLabelInterpreter',interpreter,'LineWidth',linewidth);
else
	set(axishandle,'Box','On','LineWidth',2,'FontName',fontname,...
		'FontSize',fontsize,'TickLabelInterpreter',interpreter,...
		'LineWidth',linewidth);
end
h = findobj(axishandle,'Type','line');
set(h,'LineWidth',linewidth);
set(axishandle,'Layer','top')
% if(~isempty(h))
% 	x = h.Xdata;
% 	y = h.Ydata;
% 	xlim([min(x(:)) max(x(:))]);
% 	ylim([min(y(:)) max(y(:))]);
% end

end