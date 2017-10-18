classdef customBayesianOptimization
    %BayesianOptimization    The results of a Bayesian Optimization.
    %
    %   BayesianOptimization Properties:
    %
    %     Problem definition properties:
    %       ObjectiveFcn   	- The ObjectiveFcn argument that was passed to
    %                       bayesopt.
    %       VariableDescriptions
    %                       - The VariableDescriptions argument that was
    %                       passed to bayesopt.
    %       Options         - A read-only struct containing options that
    %                       were used to perform the optimization. It has
    %                       the following fields:
    %
    %                             AcquisitionFunctionName
    %                             IsObjectiveDeterministic
    %                             ExplorationRatio
    % 
    %                             NumSeedPoints
    %                             MaxObjectiveEvaluations
    %                             MaxTime
    % 
    %                             XConstraintFcn
    %                             ConditionalVariableFcn
    %                             NumCoupledConstraints
    %                             CoupledConstraintTolerances
    %                             AreCoupledConstraintsDeterministic
    % 
    %                             Verbose
    %                             OutputFcn
    %                             SaveVariableName
    %                             SaveFileName
    %                             PlotFcn
    % 
    %                             InitialX
    %                             InitialObjective
    %                             InitialConstraintViolations
    %                             InitialErrorValues
    %                             InitialObjectiveEvaluationTimes
    %                             InitialIterationTimes
    %                             InitialUserData
    %
    %     Solution properties:
    %    	MinObjective	- The minimum observed feasible Objective
    %                       function value, with feasibility defined by the
    %                       final constraint models (including the Error
    %                       constraint model).
    %     	XAtMinObjective - A 1-by-D table representing the value of X at
    %                       the minimum observed Objective function
    %                       valuewith feasibility defined by the final
    %                       constraint models (including the Error
    %                       constraint model).
    %    	MinEstimatedObjective 
    %                       - The minimum estimated feasible Objective
    %                       function value, according to the model of
    %                       Objective, with feasibility defined by the
    %                       final constraint models (including the Error
    %                       constraint model).
    %    	XAtMinEstimatedObjective   
    %                       - A 1-by-D table representing the value of X at
    %                       the minimum estimated feasible Objective
    %                       function value, according to the model of
    %                       Objective, with feasibility defined by the
    %                       final constraint models (including the Error
    %                       constraint model).
    %       NumObjectiveEvaluations
    %                       - The number of Objective function evaluations
    %                       performed.
    %       TotalElapsedTime
    %                       - The total elapsed time of the optimization.
    %       NextPoint       - A 1-by-D table representing the next point to
    %                       be evaluated were the optimization to continue.
    %
    %     Trace properties:
    %    	XTrace          - A T-by-D table of the variable values on
    %                       which ObjectiveFcn was evaluated, where T =
    %                       NumObjectiveEvaluations and D is the number of
    %                       variables.
    %     	ObjectiveTrace  - A T-by-1 array of Objective values resulting
    %                       from the T evaluations of ObjectiveFcn on
    %                       XTrace.
    %     	ConstraintsTrace
    %                       - A T-by-K array of constraint values for K
    %                       coupled constraints resulting from the T
    %                       evaluations of ObjectiveFcn on XTrace.
    %     	UserDataTrace   - A T-by-1 cell array of the UserData outputs
    %                       of the Objective function evaluations.
    %     	ObjectiveEvaluationTimeTrace    
    %                       - A T-by-1 array of the runtime of each of
    %                    	the T evaluations of ObjectiveFcn on XTrace.
    %     	IterationTimeTrace    
    %                       - A T-by-1 array of the total time taken for
    %                       each iteration, including both function
    %                       evaluation time and overhead.
    %    	ErrorTrace      - A T-by-1 vector of error values in {-1,1}
    %                       indicating on which points in XTrace the
    %                       objective function returned a nonfinite value.
    %     	FeasibilityTrace
    %                       - A T-by-1 logical vector indicating which
    %                       points in XTrace are feasible according to the
    %                       final constraint models (including the Error
    %                       constraint model). The tolerances defined
    %                       in the 'CoupledConstraintTolerances' argument
    %                       are used to decide feasibility for each
    %                       constraint.
    %     	FeasibilityProbabilityTrace
    %                       - A T-by-1 vector of the probabilities that
    %                       each point in XTrace is feasible according to
    %                       the final constraint models (including the
    %                       Error constraint model).
    %     	IndexOfMinimumTrace    
    %                       - A T-by-1 array of integer indices indicating
    %                       which element of the trace had the minimum
    %                       feasible Objective up to that iteration, with
    %                       feasibility defined by the Constraint models at
    %                       that iteration (including the Error constraint
    %                       model).
    %     	ObjectiveMinimumTrace       
    %                       - A T-by-1 array of the minimum feasible
    %                       Objective value found up to that iteration,
    %                       with feasibility defined by the Constraint
    %                       models at that iteration (including the Error
    %                       constraint model).
    %     	EstimatedObjectiveMinimumTrace
    %                       - A T-by-1 array of the estimated minimum
    %                       Objective value found up to that iteration,
    %                       with feasibility defined by the Constraint
    %                       models at that iteration (including the Error
    %                       constraint model).
    %
    %   BayesianOptimization methods:
    %       bestPoint       - Return the best point according to a 
    %                       criterion.
    %       plot            - Create plots from a set of plot functions.
    %       predictObjective - Predict the Objective at a set of points.
    %       predictObjectiveEvaluationTime - Predict the Objective function
    %                       evaluation time at a set of points.
    %       predictConstraints - Predict Constraint values at a set of
    %                       points.
    %       predictError    - Predict the Error value at a set of points.
    %       resume          - Resume a Bayesian optimization with modified
    %                       settings and stopping criteria.
    %
    %   See also: BAYESOPT, OPTIMIZABLEVARIABLE
    
    %   Copyright 2016 The MathWorks, Inc.
    
    %% User API
    properties(Dependent, SetAccess=protected)
        % Problem definition
        ObjectiveFcn;
        VariableDescriptions;
        Options;
        % Solution / State of optimization
        MinObjective;
        XAtMinObjective;
        MinEstimatedObjective;
        XAtMinEstimatedObjective;
        NumObjectiveEvaluations;
        TotalElapsedTime;
        NextPoint;
        % Traces
        XTrace;
        ObjectiveTrace;
        ConstraintsTrace;
        UserDataTrace;
        ObjectiveEvaluationTimeTrace;
        IterationTimeTrace;
        ErrorTrace;
        FeasibilityTrace;
        FeasibilityProbabilityTrace;
        IndexOfMinimumTrace;
        ObjectiveMinimumTrace;
        EstimatedObjectiveMinimumTrace;
    end
    
    properties
        customOpts;
    end
    
    methods
        function [X, CriterionValue, Iteration] = bestPoint(this, varargin)
            % BESTPOINT Return the best point found in a Bayesian
            % optimization according to a criterion.
            %
            %   (1) X = bestPoint(RESULTS) returns the 1-row table X
            %       representing the best feasible point according to the
            %       default criterion.
            %
            %   (2) [X, CriterionValue] = bestPoint(RESULTS) also returns
            %       the value of the specified criterion at X.
            %
            %   (3) [X, CriterionValue, Iteration] = bestPoint(RESULTS)
            %       when Criterion is 'min-observed', 'min-visited-mean',
            %       or 'min-visited-upper-confidence-interval', also
            %       returns the iteration number at which the best point
            %       occurred.
            %
            %   (4) [...] = bestPoint(RESULTS, Name, Value,...)
            %       specifies additional name/value parameters:
            %
            %       Criterion   - Specifies the criterion used to choose
            %                   the best point. Possible values are the
            %                   following (case insensitive, hyphens are
            %                   optional, and unique prefixes are
            %                   recognized):
            %           'min-observed'  - X is the feasible point at which
            %                           the minimum observed Objective was
            %                           found.
            %           'min-mean'      - X is the feasible point at which
            %                           the mean of the Objective function
            %                           model is minimized.
            %           'min-upper-confidence-interval'
            %                           - X is the feasible point that
            %                           minimizes an upper confidence
            %                           interval according to the Objective
            %                           function model. See the ALPHA
            %                           parameter.
            %           'min-visited-mean'
            %                           - X is the feasible point at which
            %                           the mean of the Objective function
            %                           model is minimized, among those
            %                           points at which the Objective
            %                           function was evaluated.
            %           'min-visited-upper-confidence-interval'
            %                           - X is the feasible point that
            %                           minimizes an upper confidence
            %                           interval according to the Objective
            %                           function model, among those points
            %                           at which the Objective function was
            %                           evaluated. See the ALPHA parameter.
            %           Default: 'min-visited-upper-confidence-interval'
            %
            %       Alpha   - A probability value defining the upper
            %               confidence interval when CRITERION is
            %               'min-upper-confidence-interval' or
            %               'min-visited-upper-confidence-interval'. The
            %               upper confidence interval at X is the value Y
            %               such that Prob(mean(ObjectiveFcn(X))>Y)=Alpha
            %               according to the Objective function model.
            %               Default: .01
            
            % Parse Criterion and Alpha
            [Criterion, Alpha] = internal.stats.parseArgs({'Criterion', 'Alpha'}, ...
                {'minvisitedupperconfidenceinterval', .01}, varargin{:});
            % Parse Criterion value
            Crit = bayesoptim.parseArgValue(Criterion, {'minobserved', ...
                'minmean', 'minupperconfidenceinterval', 'minvisitedmean',...
                'minvisitedupperconfidenceinterval'});
            % Check Alpha
            if ~(isscalar(Alpha) && isnumeric(Alpha) && isreal(Alpha) && Alpha>0 && Alpha<1)
                bayesoptim.err('BadAlpha');
            end
            % Apply criterion
            Iteration = NaN;
            switch Crit
                case {'minobserved'}
                    [X, CriterionValue, Iteration] = minObservedPoint(this);
                case {'minmean'}
                    [X, CriterionValue] = minMeanPoint(this);
                case {'minvisitedmean'}
                    [X, CriterionValue, Iteration] = minVisitedMeanPoint(this);
                case {'minupperconfidenceinterval'}
                    [X, CriterionValue] = minUCIPoint(this, Alpha);
                case {'minvisitedupperconfidenceinterval'}
                    [X, CriterionValue, Iteration] = minVisitedUCIPoint(this, Alpha);
            end
        end
            
        function plot(this, varargin)
            %PLOT Plot Bayesian Optimization results
            %   (1) plot(RESULTS) or plot(RESULTS, 'all') calls all
            %       pre-defined plot functions on RESULTS.
            %
            %   (2) plot(RESULTS, plotFcn1, plotFcn2, ...) calls plot
            %       functions plotFcn1, plotFcn2, ... on RESULTS.
            %
            % Pre-defined plot functions:
            %     @plotObjectiveModel
            %     @plotAcquisitionFunction
            %     @plotObjectiveEvaluationTimeModel
            %     @plotConstraintModels
            %     @plotObjective
            %     @plotObjectiveEvaluationTime
            %     @plotMinObjective
            %     @plotElapsedTime
            if isempty(varargin) || isequal(varargin{1}, 'all')
                PlotFcn = {...
                    @plotObjectiveModel,...
                    @plotAcquisitionFunction,...
                    @plotObjectiveEvaluationTimeModel,...
                    @plotConstraintModels,...
                    @plotObjective,...
                    @plotObjectiveEvaluationTime,...
                    @plotMinObjective,...
                    @plotElapsedTime};
            else
                PlotFcn = varargin;
            end
            if ~all(cellfun(@(x)isa(x, 'function_handle'), PlotFcn))
                bayesoptim.err('PlotArgs');
            end
            for i = 1:numel(PlotFcn)
                PlotFcn{i}(this, 'standalone');
            end
            drawnow;
        end
        
        function [ConstraintViolations, SDs] = predictConstraints(this, XTable)
            %PREDICTCONSTRAINTS Predict Constraint Violations at a
            %set of points.
            %   [ConstraintViolations, SDs] =
            %   predictConstraints(Results, XTable) returns the
            %   posterior means and standard deviations of all Constraints
            %   at the points in XTable, according to a Gaussian Process
            %   model of each constraint. ConstraintViolations and SDs are
            %   both N-by-K matrices, given N rows of XTable and K coupled
            %   constraints.
            XTable = checkAndPrepareTableForPrediction(this, XTable, 'predictConstraints');
            if isempty(XTable) || this.PrivOptions.NumCoupledConstraints == 0
                ConstraintViolations = [];
                SDs = [];
            else
                for k = this.PrivOptions.NumCoupledConstraints : -1 : 1
                    if isempty(this.ConstraintGPs{k})
                        ConstraintViolations(:,k) = NaN(height(XTable),1);
                        SDs(:,k) = NaN(height(XTable),1);
                    else
                        [ConstraintViolations(:,k), SDs(:,k)] = predict(...
                            this.ConstraintGPs{k}, transformPoints(this, XTable));
                    end
                end
            end
        end
        
        function [Error, SD] = predictError(this, XTable)
            %PREDICTERROR Predict Error value at a set of points.
            %   [Error, SD] = predictError(Results, XTable) returns
            %   the posterior mean and standard deviation of the Error
            %   value at the points in XTable, according to a Gaussian
            %   Process model. Error and SD are both N-by-1 vectors, given
            %   N rows of XTable.
            XTable = checkAndPrepareTableForPrediction(this, XTable, 'predictError');
            if isempty(XTable)
                Error = [];
                SD = [];
            elseif all(this.ErrorTrace < 0)
                Error = -ones(height(XTable),1);
                SD = zeros(height(XTable),1);
            elseif isempty(this.ErrorGP)
                Error = NaN(height(XTable),1);
                SD = NaN(height(XTable),1);
            else
                [Error, SD] = predict(this.ErrorGP, transformPoints(this, XTable));
            end
        end
        
        function [Objective, SD] = predictObjective(this, XTable)
            %PREDICTOBJECTIVE Predict the Objective function at a
            %set of points.
            %   [Objective, SD] = predictObjective(Results, XTable)
            %   returns the posterior mean and standard deviation of the
            %   objective function at the points in XTable, according to a
            %   Gaussian Process model. Objective and SD are both N-by-1
            %   vectors, given N rows of XTable.
            XTable = checkAndPrepareTableForPrediction(this, XTable, 'predictObjective');
            if isempty(XTable)
                Objective = [];
                SD = [];
            elseif isempty(this.ObjectiveFcnGP)
                Objective = NaN(height(XTable),1);
                SD = NaN(height(XTable),1);
            else
                [Objective, SD] = predict(this.ObjectiveFcnGP, transformPoints(this, XTable));
            end
        end
        
        function ObjectiveEvaluationTime = predictObjectiveEvaluationTime(this, XTable)
            %PREDICTOBJECTIVEEVALUATIONTIME Predict the objective
            %function evaluation runtime at a set of points.
            %   ObjectiveEvaluationTime =
            %   predictObjectiveEvaluationTime(Results, XTable)
            %   returns the posterior mean of ObjectiveEvaluationTime at
            %   the points in XTable, according to a Gaussian Process
            %   model. ObjectiveEvaluationTime is an N-by-1 vector, given N
            %   rows of XTable.
            XTable = checkAndPrepareTableForPrediction(this, XTable, 'predictObjectiveEvaluationTime');
            if isempty(XTable)
                ObjectiveEvaluationTime = [];
            elseif isempty(this.ObjectiveEvaluationTimeGP)
                ObjectiveEvaluationTime = NaN(height(XTable),1);
            else
                logObjectiveEvaluationTime = predict(this.ObjectiveEvaluationTimeGP, transformPoints(this, XTable));
                ObjectiveEvaluationTime = exp(logObjectiveEvaluationTime);
            end
        end
        
        function NewResults = resume(Results, varargin)
            %RESUME Resume a Bayesian Optimization.
            %   RESULTS = resume(RESULTS, 'PARAM1', val1, 'PARAM2'
            %   ,val2,...) continues running an optimization with modified
            %   options until new stopping criteria are met. All name/value
            %   parameters accepted by bayesopt are accepted, except
            %   'NumSeedPoints' and the 'Initial___' parameters. Important
            %   parameters for resuming are:
            %       'MaxObjectiveEvaluations' 
            %                   - Specifies the additional number of
            %                     function evaluations to run. Default: 30
            %       'MaxTime'   - Specifies the additional number of
            %                     seconds to run. Default: Inf
            %       'VariableDescriptions'
            %                   - Specifies new variable descriptions to
            %                     use. The variable names to be optimized
            %                     must remain the same. Numeric variables
            %                     may have their Range, Type, and Transform
            %                     changed, but must remain numeric.
            %                     Categorical variables may not be changed.
            if iAnyInitializationArgs(varargin)
                bayesoptim.err('NoInitIfResume');
            else
                % Parse out new variable descriptions
                [PassedVariableDescrips, RemainingArgs] = iParseVariableDescriptionsArg(varargin);
                if isempty(PassedVariableDescrips)
                    VariableDescrips = Results.PrivOptions.VariableDescriptions;
                else
                    checkVariableDescriptions(PassedVariableDescrips, Results.PrivOptions.VariableDescriptions);
                    VariableDescrips = PassedVariableDescrips;
                end
                % Create NewOptions, install new variable descriptions
                NewOptions = Results.PrivOptions;
                NewOptions.VariableDescriptions = VariableDescrips;
                % Create an options obj just to get new stopping criteria
                ParsedOptions = bayesoptim.BayesoptOptions(Results.PrivOptions.ObjectiveFcn, VariableDescrips, RemainingArgs);
                % Interpret new stopping criteria as increments to the current state
                NewOptions.MaxObjectiveEvaluations = Results.NumObjectiveEvaluations + ParsedOptions.MaxObjectiveEvaluations;
                NewOptions.MaxTime = Results.TotalElapsedTime + ParsedOptions.MaxTime;
                % Save current XTrace
                CurXTrace = Results.XTrace;
                % Install new options (this would change XTrace)
                Results.PrivOptions = NewOptions;
                % Update XTrain because variableDescriptions may have changed
                Results.XTrain = XTraceToXTrain(Results, CurXTrace);
                % Run the optimization using the new options
                NewResults = run(Results);
            end
        end
    end
    
    %% Hidden API
    methods(Hidden)
        function this = customBayesianOptimization(Options, customOpts)
            this.PrivOptions = Options;
            this.customOpts = customOpts;
            this = run(this);
        end
        
        % Output functions
        function Stop = assignInBase(bo, State)
            %assignInBase     OutputFcn for assigning a Bayesian
            %Optimization to a variable in the base workspace. 
            %   Stop = assignInBase(bo, State) assigns the
            %   BayesianOptimization instance 'bo' to a variable in the
            %   base workspace. The variable name is obtained from
            %   bo.Options.SaveVariableName.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            narginchk(2,2);
            Stop = false;
            try
                assignin('base', bo.PrivOptions.SaveVariableName, bo);
            catch me
                bayesoptim.warn('AssignInBaseFailed');
                disp(me);
            end
        end
        
        function Stop = saveToFile(bo, State)
            %saveToFile     OutputFcn for saving a Bayesian Optimization to
            %a file. 
            %   Stop = saveToFile(bo, State) saves the BayesianOptimization
            %   instance 'bo' to a file. The file name is obtained from
            %   bo.Options.SaveFileName.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            narginchk(2,2);
            Stop = false;
            SaveFilename = bo.PrivOptions.SaveFileName;
            prevFilename = sprintf('%s.PREV', SaveFilename);
            try
                if exist(['./' SaveFilename], 'file')
                    if exist(['./' prevFilename], 'file')
                        delete(prevFilename);
                    end
                    copyfile(SaveFilename, prevFilename);
                end
                BayesoptResults = bo;
                save(SaveFilename, 'BayesoptResults');
            catch me
                bayesoptim.warn('SaveToFileFailed', SaveFilename);
                disp(me);
            end
        end
        
        % Predefined plot functions
        function Stop = plotAcquisitionFunction(bo, State)
            %plotAcquisitionFunction     PlotFcn for plotting the
            %Acquisition Function surface.
            %   Stop = plotAcquisitionFunction(bo, State) plots the
            %   Acquisition function surface. The value is set to NaN for
            %   points that violate the xConstraintFcn.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curves;
            narginchk(2,2);
            Stop = false;
            if bo.NumVars > 2 || ismember(bo.PrivOptions.AcquisitionFunctionName, {'grid','random'})
                return;
            end
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    Axes = setupPlot(bo);
                    Curves = {};
                case {'iteration', 'done'}
                    if isvalid(Axes)
                        [Axes, Curves] = updatePlot(bo, Axes, Curves);
                    end
                case 'standalone'
                    StandaloneAxes = setupPlot(bo);
                    updatePlot(bo, StandaloneAxes, {});
            end
            
            function Axes = setupPlot(Results)
                screenSize = get(groot,'ScreenSize');
                L = screenSize(1);
                B = screenSize(2);
                W = screenSize(3);
                H = screenSize(4);
                left = W/3 + 1;
                bottom = H/2 - 100;
                width = W/3 - 50;
                height = H/2;
                f = figure('Position',[left, bottom, width, height],'Tag','bayesopt.AcqFcn');
                Axes = axes(f);
                title(Axes, bayesoptim.infoString('AcqFcn', Results.PrivOptions.AcquisitionFunctionName));
                if Results.NumVars == 2
                    view(Axes, -45, 15);
                end
            end
            
            function [Axes, Curves] = updatePlot(Results, Axes, Curves)
                import bayesoptim.*
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                FModel = Results.ObjectiveFcnGP;
                if ismember(Options.AcquisitionFunctionName, {'grid', 'random'})
                    return;
                end
                if isempty(FModel)
                    return;
                end
                ObjectiveEvaluationTimeModel = Results.ObjectiveEvaluationTimeGP;
                if isempty(ObjectiveEvaluationTimeModel) && ismember(Options.AcquisitionFunctionName, ...
                        {'expectedimprovementpersecond', 'expectedimprovementpersecondplus'})
                    return;
                end
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = table2array(untransformPoints(Results, EvalGrid, false));
                        % Compute AF at gridpoints and at next point
                        switch Options.AcquisitionFunctionName
                            case 'probabilityofimprovement'
                                AFGrid = probabilityOfImprovement(LegalizedEvalGrid, ...
                                    FModel, Results.IncumbentF, Results.ObjectiveFcnGP.Sigma);
                                AFNext = probabilityOfImprovement(legalizePoints(Results, Results.XNext), ...
                                    FModel, Results.IncumbentF, Results.ObjectiveFcnGP.Sigma);
                            case 'lowerconfidencebound'
                                AFGrid = lowerConfidenceBound(LegalizedEvalGrid, ...
                                    FModel, 2);
                                AFNext = lowerConfidenceBound(legalizePoints(Results, Results.XNext), ...
                                    FModel, 2);
                            case {'expectedimprovement', 'expectedimprovementplus'}
                                AFGrid = expectedImprovement(LegalizedEvalGrid, ...
                                    FModel, ...
                                    Results.IncumbentF);
                                AFNext = expectedImprovement(legalizePoints(Results, Results.XNext), ...
                                    FModel, ...
                                    Results.IncumbentF);
                            case {'expectedimprovementpersecond', 'expectedimprovementpersecondplus'}
                                AFGrid = expectedImprovementPerSecond(LegalizedEvalGrid, ...
                                    FModel, ObjectiveEvaluationTimeModel, Results.IncumbentF);
                                AFNext = expectedImprovementPerSecond(legalizePoints(Results, Results.XNext), ...
                                    FModel, ObjectiveEvaluationTimeModel, Results.IncumbentF);
                        end
                        % NaN-out points that violate xConstraint
                        mask = double(satisfiesXConstraint(Results, EvalGrid));
                        mask(mask==0) = NaN;
                        AFGrid = AFGrid .* mask;
                        % Apply constraint weighting
                        AFGrid = AFGrid .* ProbAllConstraintsSatisfied(Results, EvalGrid);
                        % Plot surface
                        Curves{end+1} = plot(Axes, PlotGrid, AFGrid, 'r');
                        % Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot(Axes, XNextToPlot, AFNext, MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, Options.AcquisitionFunctionName);
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        C = bayesoptim.suppressWarnings();
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = table2array(untransformPoints(Results, EvalGrid, false));
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % Compute AF at gridpoints and at next point
                        switch Options.AcquisitionFunctionName
                            case 'probabilityofimprovement'
                                AFGrid = probabilityOfImprovement(LegalizedEvalGrid, ...
                                    FModel, Results.IncumbentF, Results.ObjectiveFcnGP.Sigma);
                                AFNext = probabilityOfImprovement(legalizePoints(Results, Results.XNext), ...
                                    FModel, Results.IncumbentF, Results.ObjectiveFcnGP.Sigma);
                            case 'lowerconfidencebound'
                                AFGrid = lowerConfidenceBound(LegalizedEvalGrid, ...
                                    FModel, 2);
                                AFNext = lowerConfidenceBound(legalizePoints(Results, Results.XNext), ...
                                    FModel, 2);
                            case {'expectedimprovement', 'expectedimprovementplus'}
                                AFGrid = expectedImprovement(LegalizedEvalGrid, ...
                                    FModel, Results.IncumbentF);
                                AFNext = expectedImprovement(legalizePoints(Results, Results.XNext), ...
                                    FModel, Results.IncumbentF);
                            case {'expectedimprovementpersecond', 'expectedimprovementpersecondplus'}
                                AFGrid = expectedImprovementPerSecond(LegalizedEvalGrid, ...
                                    FModel, ObjectiveEvaluationTimeModel, Results.IncumbentF);
                                AFNext = expectedImprovementPerSecond(legalizePoints(Results, Results.XNext), ...
                                    FModel, ObjectiveEvaluationTimeModel, Results.IncumbentF);
                        end
                        % NaN-out points that violate xConstraint
                        mask = double(satisfiesXConstraint(Results, LegalizedEvalGrid));
                        mask(mask==0) = NaN;
                        AFGrid = AFGrid .* mask;
                        % Apply constraint weighting
                        AFGrid = AFGrid .* ProbAllConstraintsSatisfied(Results, LegalizedEvalGrid);
                        
                        % Plot surface
                        Curves{end+1} = surfc(Axes, X1Plot, X2Plot, reshape(AFGrid, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'b', 'FaceAlpha', .6);
                        % Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), AFNext, ...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        zlabel(Axes, Options.AcquisitionFunctionName);
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                end
            end
        end
        
        function Stop = plotConstraintModels(bo, State)
            %plotConstraintModels     PlotFcn for plotting Constraint-model
            %surfaces. 
            %   Stop = plotConstraintModels(bo, State) plots each
            %   Constraint-model surface, including the built-in Error
            %   constraint if there were any errors. It also plots a
            %   Prob(Feasible) surface.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curves NoErrorsYet;
            narginchk(2,2);
            Stop = false;
            if ~bo.PrivOptions.FitModels || bo.NumVars > 2
                return;
            end
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            NumConstraints = bo.PrivOptions.NumCoupledConstraints;
            % Order of plots is: {C1, ..., Cn, ProbFeas, Error}. Note that
            % the error plot will not exist until there is an error.
            switch State
                case 'initial'
                    % User constraints
                    for i = 1:NumConstraints
                        Axes{i} = setupPlot(bo, bayesoptim.infoString('ConstraintNum', i), ...
                            bayesoptim.infoString('DegreeOfViolation'), ...
                            sprintf('bayesopt.Constraint%d', i));
                        Curves{i} = {};
                    end
                    % Prob feasible
                    Axes{NumConstraints+1} = setupPlot(bo, bayesoptim.infoString('ProbFeas'), ...
                        bayesoptim.infoString('ProbFeas'), 'bayesopt.ProbFeas');
                    Curves{NumConstraints+1} = {};
                    % Error Constraint
                    NoErrorsYet = true;
                    if ~isempty(bo.ErrorGP)
                        NoErrorsYet = false;
                        Axes{NumConstraints+2} = setupPlot(bo, bayesoptim.infoString('ErrorConstraint'), ...
                            bayesoptim.infoString('DegreeOfViolation'), 'bayesopt.ErrorConstraint');
                        Curves{NumConstraints+2} = {};
                    end
                case {'iteration', 'done'}
                    % User constraints
                    for i = 1:NumConstraints
                        if isvalid(Axes{i})
                            [Axes{i}, Curves{i}] = updateConstraintPlot(bo, i, Axes{i}, Curves{i});
                        end
                    end
                    % Prob feasible
                    if isvalid(Axes{NumConstraints+1})
                        [Axes{NumConstraints+1}, Curves{NumConstraints+1}] = updateProbFeasPlot(bo, Axes{NumConstraints+1}, Curves{NumConstraints+1});
                    end
                    % Error Constraint
                    if ~isempty(bo.ErrorGP) 
                        if NoErrorsYet
                            NoErrorsYet = false;
                            Axes{NumConstraints+2} = setupPlot(bo, bayesoptim.infoString('ErrorConstraint'), ...
                                bayesoptim.infoString('DegreeOfViolation'), 'bayesopt.ErrorConstraint');
                            Curves{NumConstraints+2} = {};
                        end
                        if isvalid(Axes{NumConstraints+2})
                            [Axes{NumConstraints+2}, Curves{NumConstraints+2}] = updateErrorPlot(...
                                bo, Axes{NumConstraints+2}, Curves{NumConstraints+2});
                        end
                    end
                case 'standalone'
                    % User constraints
                    for i = 1:NumConstraints
                        StandaloneAxes{i} = setupPlot(bo, bayesoptim.infoString('ConstraintNum', i), ...
                            bayesoptim.infoString('DegreeOfViolation'), ...
                            sprintf('bayesopt.Constraint%d', i));
                        updateConstraintPlot(bo, i, StandaloneAxes{i}, {});
                    end
                    % Prob feasible
                    StandaloneAxes{NumConstraints+1} = setupPlot(bo, bayesoptim.infoString('ProbFeas'), ...
                        bayesoptim.infoString('ProbFeas'), 'bayesopt.ProbFeas');
                    updateProbFeasPlot(bo, StandaloneAxes{NumConstraints+1}, {});
                    % Error Constraint
                    if ~isempty(bo.ErrorGP)
                        StandaloneAxes{NumConstraints+2} = setupPlot(bo, bayesoptim.infoString('ErrorConstraint'), ...
                            bayesoptim.infoString('DegreeOfViolation'), 'bayesopt.ErrorConstraint');
                        updateErrorPlot(bo, StandaloneAxes{NumConstraints+2}, {});
                    end
            end
            
            function Axes = setupPlot(Results, Title, ZLabel, Tag)
                screenSize = get(groot,'ScreenSize');
                L = screenSize(1);
                B = screenSize(2);
                W = screenSize(3);
                H = screenSize(4);
                left = L;
                bottom = H/2 - 100;
                width = W/3 - 50;
                height = H/2;
                f = figure('Position',[left, bottom, width, height],'Tag',Tag);
                Axes = axes(f);
                title(Axes, Title);
                zlabel(Axes, ZLabel)
                if Results.NumVars == 2
                    view(Axes, -45, 15);
                end
            end
            
            function [Axes, Curves] = updateConstraintPlot(Results, ConstraintNum, Axes, Curves)
                import bayesoptim.*
                if isempty(Results.ConstraintGPs{ConstraintNum})
                    return;
                end
                
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                Model = Results.ConstraintGPs{ConstraintNum};
                XTrace = Results.XTrace;
                ValTrace = Results.ConstraintsTrace(:,ConstraintNum);
                
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,1};
                        % (1) Eval GPRs on eval grid
                        YPred = predict(Model, LegalizedEvalGrid);
                        % (2) Plot training points on all axes
                        plot(Axes, XPlottable, ValTrace, 'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot Model
                        Curves{end+1} = plot(Axes, PlotGrid, YPred, 'r');
                        Curves{end+1} = plot(Axes, [PlotGrid(1); PlotGrid(end)], [0;0], 'g');     % Plot line at Y=0
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        ConstraintVals = predictConstraints(Results, Results.NextPoint);
                        Curves{end+1} = plot(Axes, XNextToPlot, ConstraintVals(ConstraintNum), ...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, bayesoptim.infoString('ConstraintNum', ConstraintNum));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        C = bayesoptim.suppressWarnings();
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,:};
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % (1) Eval GPRs on eval grid
                        YPred = predict(Model, LegalizedEvalGrid);
                        % (2) Plot training points
                        Curves{end+1} = plot3(Axes, XPlottable(:,1), XPlottable(:,2), ValTrace, ...
                            'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot model
                        Curves{end+1} = surfc(Axes, X1Plot, X2Plot, reshape(YPred, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'r', 'FaceAlpha', .5, 'LineStyle', '-');
                        Curves{end+1} = surf(Axes, X1Plot, X2Plot, zeros(reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'g', 'FaceAlpha', .5, 'LineStyle','-');     % Plot surface at Y=0
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        ConstraintVals = predictConstraints(Results, Results.NextPoint);
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), ConstraintVals(ConstraintNum),...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        zlabel(Axes, bayesoptim.infoString('DegreeOfViolation'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                        
                end
            end
            
            function [Axes, Curves] = updateErrorPlot(Results, Axes, Curves)
                import bayesoptim.*
                if isempty(Results.ErrorGP)
                    return;
                end
                
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                Model = Results.ErrorGP;
                XTrace = Results.XTrace;
                ValTrace = Results.ErrorTrace;
                
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,1};
                        % (1) Eval GPRs on eval grid
                        YPred = predict(Model, LegalizedEvalGrid);
                        % (2) Plot training points
                        plot(Axes, XPlottable, ValTrace, 'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot Model
                        Curves{end+1} = plot(Axes, PlotGrid, YPred, 'r');
                        Curves{end+1} = plot(Axes, [PlotGrid(1); PlotGrid(end)], [0;0], 'g');     % Plot line at Y=0
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot(Axes, XNextToPlot, predictError(Results, Results.NextPoint), ...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, bayesoptim.infoString('DegreeOfViolation'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        C = bayesoptim.suppressWarnings();
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,:};
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % (1) Eval GPRs on eval grid
                        YPred = predict(Model, LegalizedEvalGrid);
                        % (2) Plot training points
                        Curves{end+1} = plot3(Axes, XPlottable(:,1), XPlottable(:,2), ValTrace, ...
                            'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot model
                        Curves{end+1} = surfc(Axes, X1Plot, X2Plot, reshape(YPred, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'r', 'FaceAlpha', .5, 'LineStyle', '-');
                        Curves{end+1} = surf(Axes, X1Plot, X2Plot, zeros(reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'g', 'FaceAlpha', .5, 'LineStyle','-');     % Plot surface at Y=0
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), ...
                            predictError(Results, Results.NextPoint),...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        zlabel(Axes, bayesoptim.infoString('DegreeOfViolation'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                        
                end
            end
            
            function [Axes, Curves] = updateProbFeasPlot(Results, Axes, Curves)
                import bayesoptim.*
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                XTrace = Results.XTrace;
                ValTrace = Results.FeasibilityProbabilityTrace;
                
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,1};
                        % (1) Eval GPRs on eval grid
                        YPred = ProbAllConstraintsSatisfied(Results, LegalizedEvalGrid);
                        % (2) Plot training points on all axes
                        Curves{end+1} = plot(Axes, XPlottable, ValTrace, 'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot Model
                        Curves{end+1} = plot(Axes, PlotGrid, YPred, 'r');
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot(Axes, XNextToPlot, ...
                            ProbAllConstraintsSatisfied(Results, Results.XNext), ...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, bayesoptim.infoString('ProbFeas'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        C = bayesoptim.suppressWarnings();
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,:};
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % (1) Eval GPRs on eval grid
                        YPred = ProbAllConstraintsSatisfied(Results, LegalizedEvalGrid);
                        % (2) Plot training points
                        Curves{end+1} = plot3(Axes, XPlottable(:,1), XPlottable(:,2), ValTrace, ...
                            'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot model
                        Curves{end+1} = surfc(Axes, X1Plot, X2Plot, reshape(YPred, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'r', 'FaceAlpha', .5, 'LineStyle', '-');
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), ...
                            ProbAllConstraintsSatisfied(Results, Results.XNext),...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        zlabel(Axes, bayesoptim.infoString('ProbFeas'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                        
                end
            end
        end
        
        function Stop = plotElapsedTime(bo, State)
            %plotElapsedTime     PlotFcn for plotting total elapsed time
            %versus the number of function evaluations. 
            %   Stop = plotElapsedTime(bo, State) plots total elapsed time
            %   versus the number of function evaluations.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curves;
            narginchk(2,2);
            Title = bayesoptim.infoString('ElapsedTimeTitle');
            XLabel = bayesoptim.infoString('FEvals');
            YLabel = bayesoptim.infoString('TotElapsed');
            Stop = false;
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    f = figure('Tag','bayesopt.ElapsedTime');
                    Axes = axes(f);
                    title(Axes, Title);
                    xlabel(Axes, XLabel);
                    ylabel(Axes, YLabel);
                    Curves = {};
                case {'iteration', 'done'}
                    cellfun(@delete, Curves);
                    if isvalid(Axes)
                        obj = nancumsum(bo.ObjectiveEvaluationTimeTrace);
                        tot = nancumsum(bo.IterationTimeTrace);
                        Curves{end+1} = plotCurve(tot, Axes, 'b');
                        Curves{end+1} = plotCurve(obj, Axes, 'g');
                        Curves{end+1} = plotCurve(tot-obj, Axes, 'r');
                        legend(Axes, {bayesoptim.infoString('TotTime'),...
                            bayesoptim.infoString('FEvalTime'),...
                            bayesoptim.infoString('ModelingTime')},...
                            'Location', 'best');
                    end
                case 'standalone'
                    obj = nancumsum(bo.ObjectiveEvaluationTimeTrace);
                    tot = nancumsum(bo.IterationTimeTrace);
                    f = figure('Tag','bayesopt.ElapsedTime');
                    StandaloneAxes = axes(f);
                    title(StandaloneAxes, Title);
                    xlabel(StandaloneAxes, XLabel);
                    ylabel(StandaloneAxes, YLabel);
                    plotCurve(tot, StandaloneAxes, 'b');
                    plotCurve(obj, StandaloneAxes, 'g');
                    plotCurve(tot-obj, StandaloneAxes, 'r');
                    legend(StandaloneAxes, {bayesoptim.infoString('TotTime'),...
                        bayesoptim.infoString('FEvalTime'),...
                        bayesoptim.infoString('ModelingTime')}, ...
                        'Location', 'best');
            end
            
            function Curve = plotCurve(Trace, Axes, color)
                HoldOff = holdOn(Axes);
                Curve = plot(Axes, 1:numel(Trace), Trace, ['-' color 'o'], 'MarkerSize', 2);
            end
        end
        
        function Stop = plotMinObjective(bo, State)
            %plotMinObjective     PlotFcn for plotting the minimum observed
            %Objective versus the number of function evaluations. 
            %   Stop = plotMinObjective(bo, State) plots the minimum
            %   observed Objective versus the number of function
            %   evaluations.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes1 Curves1;
            narginchk(2,2);
            Title1 = bayesoptim.infoString('MinObjTitle');
            XLabel1 = bayesoptim.infoString('FEvals');
            YLabel = bayesoptim.infoString('MinObj');
            Stop = false;
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    f = figure('Tag','bayesopt.MinObjective');
                    Axes1 = axes(f);
                    title(Axes1, Title1);
                    xlabel(Axes1, XLabel1);
                    ylabel(Axes1, YLabel);
                    Axes1.YAxisLocation='right';
                    Curves1 = {};
                case {'iteration', 'done'}
                    cellfun(@delete, Curves1);
                    if isvalid(Axes1)
                        Curves1{end+1} = plotCurve(1:numel(bo.ObjectiveMinimumTrace), ...
                            bo.ObjectiveMinimumTrace, Axes1, 'b');
                        if all(isnan(bo.EstimatedObjectiveMinimumTrace))
                            legend(Axes1, {bayesoptim.infoString('MinObsObj')}, 'Location', 'best');
                        else
                            Curves1{end+1} = plotCurve(1:numel(bo.ObjectiveMinimumTrace), ...
                                bo.EstimatedObjectiveMinimumTrace, Axes1, 'g');
                            legend(Axes1, {bayesoptim.infoString('MinObsObj'), ...
                                bayesoptim.infoString('EstMinObj')}, 'Location', 'best');
                        end
                    end
                case 'standalone'
                    f = figure('Tag','bayesopt.MinObjective');
                    SAxes1 = axes(f);
                    title(SAxes1, Title1);
                    xlabel(SAxes1, XLabel1);
                    ylabel(SAxes1, YLabel);
                    SAxes1.YAxisLocation='right';
                    plotCurve(1:numel(bo.ObjectiveMinimumTrace), bo.ObjectiveMinimumTrace, SAxes1, 'b');
                    if all(isnan(bo.EstimatedObjectiveMinimumTrace))
                        legend(SAxes1, {bayesoptim.infoString('MinObsObj')}, 'Location', 'best');
                    else
                        plotCurve(1:numel(bo.ObjectiveMinimumTrace), bo.EstimatedObjectiveMinimumTrace, SAxes1, 'g');
                        legend(SAxes1, {bayesoptim.infoString('MinObsObj'), bayesoptim.infoString('EstMinObj')},...
                            'Location', 'best');
                    end
            end
            
            function Curve = plotCurve(X, TraceToPlot, Axes, color)
                HoldOff = holdOn(Axes);
                Curve = plot(Axes, X, TraceToPlot, ['-' color 'o'], 'MarkerSize', 2);
            end
        end
        
        function Stop = plotObjective(bo, State)
            %plotObjective     PlotFcn for plotting each observed Objective
            %function value versus the number of function evaluations. 
            %   Stop = plotObjective(bo, State) plots each observed
            %   Objective function value versus the number of function
            %   evaluations.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curve;
            narginchk(2,2);
            Title = bayesoptim.infoString('ObjFunTitle');
            XLabel = bayesoptim.infoString('FEvals');
            YLabel = bayesoptim.infoString('ObjFun');
            Stop = false;
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    f = figure('Tag','bayesopt.Objective');
                    Axes = axes(f);
                    title(Axes, Title);
                    xlabel(Axes, XLabel);
                    ylabel(Axes, YLabel);
                    Curve = [];
                case {'iteration', 'done'}
                    delete(Curve);
                    if isvalid(Axes)
                        Curve = plotCurve(bo.ObjectiveTrace, Axes);
                    end
                case 'standalone'
                    f = figure('Tag','bayesopt.Objective');
                    StandaloneAxes = axes(f);
                    title(StandaloneAxes, Title);
                    xlabel(StandaloneAxes, XLabel);
                    ylabel(StandaloneAxes, YLabel);
                    plotCurve(bo.ObjectiveTrace, StandaloneAxes);
            end
            
            function Curve = plotCurve(ObjectiveTrace, Axes)
                HoldOff = holdOn(Axes);
                Curve = plot(Axes, 1:numel(ObjectiveTrace), ObjectiveTrace, '-bo', 'MarkerSize', 2);
            end
        end
        
        function Stop = plotObjectiveEvaluationTime(bo, State)
            %plotObjectiveEvaluationTime     PlotFcn for plotting each
            %observed Objective fucntion evaluation time versus the number
            %of function evaluations. 
            %   Stop = plotObjectiveEvaluationTime(bo, State) plots each
            %   observed Objective function evaluation time versus the
            %   number of function evaluations.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curve;
            narginchk(2,2);
            Title = bayesoptim.infoString('FEvalTimeTitle');
            XLabel = bayesoptim.infoString('FEvals');
            YLabel = bayesoptim.infoString('FEvalTime');
            Stop = false;
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    f = figure('Tag','bayesopt.ObjectiveEvaluationTime');
                    Axes = axes(f);
                    title(Axes, Title);
                    xlabel(Axes, XLabel);
                    ylabel(Axes, YLabel);
                    Curve = [];
                case {'iteration', 'done'}
                    delete(Curve);
                    if isvalid(Axes)
                        Curve = plotCurve(bo.ObjectiveEvaluationTimeTrace, Axes);
                    end
                case 'standalone'
                    f = figure('Tag','bayesopt.ObjectiveEvaluationTime');
                    StandaloneAxes = axes(f);
                    title(StandaloneAxes, Title);
                    xlabel(StandaloneAxes, XLabel);
                    ylabel(StandaloneAxes, YLabel);
                    plotCurve(bo.ObjectiveEvaluationTimeTrace, StandaloneAxes);
            end
            
            function Curve = plotCurve(Trace, Axes)
                HoldOff = holdOn(Axes);
                Curve = plot(Axes, 1:numel(Trace), Trace, '-bo', 'MarkerSize', 2);
            end
        end

        function Stop = plotObjectiveEvaluationTimeModel(bo, State)
            %plotObjectiveEvaluationTimeModel     PlotFcn for plotting the
            %ObjectiveEvaluationTime-model surface. 
            %   Stop = plotObjectiveEvaluationTimeModel(bo, State) plots
            %   the ObjectiveEvaluationTime-model surface.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curves;
            narginchk(2,2);
            Stop = false;
            if ~bo.PrivOptions.FitModels || bo.NumVars > 2
                return;
            end
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    Axes = setupPlot(bo);
                    Curves = {};
                case {'iteration', 'done'}
                    if isvalid(Axes)
                        [Axes, Curves] = updatePlot(bo, Axes, Curves);
                    end
                case 'standalone'
                    StandaloneAxes = setupPlot(bo);
                    updatePlot(bo, StandaloneAxes, {});
            end
            
            function Axes = setupPlot(Results)
                screenSize = get(groot,'ScreenSize');
                L = screenSize(1);
                B = screenSize(2);
                W = screenSize(3);
                H = screenSize(4);
                left = L;
                bottom = H/2 - 100;
                width = W/3 - 50;
                height = H/2;
                f = figure('Position',[left, bottom, width, height], 'Tag','bayesopt.ObjectiveEvaluationTimeModel');
                Axes = axes(f);
                title(Axes, bayesoptim.infoString('FEvalTimeModel'));
                zlabel(Axes, bayesoptim.infoString('EstFEvalTime'));
                if Results.NumVars == 2
                    view(Axes, -45, 15);
                end
            end
            
            function [Axes, Curves] = updatePlot(Results, Axes, Curves)
                import bayesoptim.*
                if isempty(Results.ObjectiveEvaluationTimeGP)
                    return;
                end
                
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                Model = Results.ObjectiveEvaluationTimeGP;
                XTrace = Results.XTrace;
                ValTrace = Results.ObjectiveEvaluationTimeTrace;
                ValTrace(isnan(Results.ObjectiveTrace)) = NaN;      % Only plot valid function evaluations
                
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,1};
                        % (1) Eval GPRs on eval grid
                        YPred = exp(predict(Model, LegalizedEvalGrid));
                        % (2) Plot training points on all axes
                        plot(Axes, XPlottable, ValTrace, 'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot Model's mean
                        Curves{end+1} = plot(Axes, PlotGrid, YPred, 'r');
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot(Axes, XNextToPlot, predictObjectiveEvaluationTime(...
                            Results, Results.NextPoint), MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, bayesoptim.infoString('FEvalTime'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        C = bayesoptim.suppressWarnings();
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,:};
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % (1) Eval GPRs on eval grid
                        YPred = exp(predict(Model, LegalizedEvalGrid));
                        % (2) Plot training points
                        Curves{end+1} = plot3(Axes, XPlottable(:,1), XPlottable(:,2), ValTrace, ...
                            'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot FModel's mean and minimum
                        Curves{end+1} = surfc(Axes, X1Plot, X2Plot, reshape(YPred, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'r', 'FaceAlpha', .5, 'LineStyle', '-');
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), ...
                            predictObjectiveEvaluationTime(Results, Results.NextPoint), ...
                            MarkerSpec{:});
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                        
                end
            end
            
        end

        function Stop = plotObjectiveModel(bo, State)
            %plotObjectiveModel     PlotFcn for plotting the Objective
            %function model surface, along with the estimated location of
            %the minimum. 
            %   Stop = plotObjectiveModel(bo, State) plots the Objective
            %   function model surface, along with the estimated location
            %   of the minimum.
            %
            % See also: BAYESOPT, BAYESIANOPTIMIZATION
            persistent Axes Curves;
            narginchk(2,2);
            Stop = false;
            VarSpec = bo.PrivOptions.VarSpec;
            if ~bo.PrivOptions.FitModels || bo.NumVars > 2
                return;
            end
            set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
            switch State
                case 'initial'
                    Axes = setupFPlot(bo);
                    Curves = {};
                case {'iteration', 'done'}
                    if isvalid(Axes)
                        [Axes, Curves] = updateFPlot(bo, Axes, Curves);
                    end
                case 'standalone'
                    StandaloneAxes = setupFPlot(bo);
                    updateFPlot(bo, StandaloneAxes, {});
            end
            
            function Axes = setupFPlot(Results)
                screenSize = get(groot,'ScreenSize');
                L = screenSize(1);
                B = screenSize(2);
                W = screenSize(3);
                H = screenSize(4);
                % (1) Plot True Objective on grid
                left = L;
                bottom = H/2 - 100;
                width = W/3 - 50;
                height = H/2;
                f = figure('Position',[left, bottom, width, height], 'Tag','bayesopt.ObjectiveModel');
                Axes = axes(f);
                title(Axes, bayesoptim.infoString('ObjModel'));
                zlabel(Axes, bayesoptim.infoString('EstObj'));
                if Results.NumVars == 2
                    view(Axes, -45, 15);
                end
            end
            
            function [Axes, Curves] = updateFPlot(Results, Axes, Curves)
                import bayesoptim.*
                set(0, 'defaultfigureunits', 'pixels');     % Make sure units are pixels
                HoldOff = holdOn(Axes);
                Options = Results.PrivOptions;
                VarSpec = Options.VarSpec;
                FModel = Results.ObjectiveFcnGP;
                if isempty(FModel)
                    return;
                end
                XTrace = Results.XTrace;
                ObjectiveTrace = Results.ObjectiveTrace;
                ModelMinXTable = untransformPoints(Results, Results.IncumbentX, true);
                ModelMinF = Results.IncumbentF;
                
                cellfun(@delete, Curves);
                Curves = {};
                switch Results.NumVars
                    case 1
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,1};
                        % (1) Eval GPRs on eval grid
                        % Objective
                        [FPred, YSD, ~] = predict(FModel, LegalizedEvalGrid);
                        FSD = funStd(YSD, FModel.Sigma);
                        % (2) Plot training points on all axes
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        Curves{end+1} = plot(Axes, XPlottable, ObjectiveTrace, 'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot FModel's mean, envelope
                        Curves{end+1} = plot(Axes, PlotGrid, FPred, 'r');
                        Curves{end+1} = plot(Axes, PlotGrid, FPred + FSD, 'b');
                        Curves{end+1} = plot(Axes, PlotGrid, FPred - FSD, 'b');
                        Curves{end+1} = plot(Axes, PlotGrid, FPred + FModel.Sigma, 'c');
                        Curves{end+1} = plot(Axes, PlotGrid, FPred - FModel.Sigma, 'c');
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot(Axes, XNextToPlot, predictObjective(Results, Results.NextPoint), ...
                            MarkerSpec{:});
                        % Plot minimum if it exists
                        if ~isempty(ModelMinXTable)
                            XMinPlottable = makeXTablePlottable(Results, ModelMinXTable);
                            Curves{end+1} = plot(Axes, XMinPlottable, ModelMinF, 'r*');
                        end
                        % Plot legend, which has either 5 or 6 items
                        LegendStrings = {...
                            bayesoptim.infoString('Observed'),...
                            bayesoptim.infoString('ModelMean'),...
                            bayesoptim.infoString('ModelErrorbars'),...
                            bayesoptim.infoString('NoiseErrorbars'),...
                            bayesoptim.infoString('NextPoint')};
                        CurveNums = [1,2,3,5,7];
                        if ~isempty(ModelMinXTable)
                            LegendStrings{end+1} = bayesoptim.infoString('ModelMin');
                            CurveNums(end+1) = 8;
                        end
                        legend(Axes, [Curves{CurveNums}], LegendStrings{:}, 'Location', 'best');
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, bayesoptim.infoString('EstObj'));
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs, VarSpec.UBs];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                    case 2
                        C = bayesoptim.suppressWarnings();
                        % Create grid to evaluate ObjectiveFcn on (in design matrix shape)
                        EvalGrid = NDGrid(VarSpec.LBTrans, VarSpec.UBTrans, Options.NumPlotGrid);
                        LegalizedEvalGrid = legalizePoints(Results, EvalGrid);
                        % Transform eval grid to native space for plotting
                        PlotGrid = untransformPoints(Results, EvalGrid, false);
                        PlotGrid = PlotGrid{:,:};
                        % Make a version of the PlotGrid in meshgrid shape:
                        X1Plot = reshape(PlotGrid(:,1), reverse(Options.NumPlotGrid));
                        X2Plot = reshape(PlotGrid(:,2), reverse(Options.NumPlotGrid));
                        % (1) Eval GPRs on eval grid
                        % Objective
                        FPred = predict(FModel, LegalizedEvalGrid);
                        % (2) Plot training points
                        XPlottable = makeXTablePlottable(Results, XTrace);
                        Curves{end+1} = plot3(Axes, XPlottable(:,1), XPlottable(:,2), ObjectiveTrace, ...
                            'bo', 'MarkerFaceColor', 'blue');
                        % (3) Plot FModel's mean
                        SurfAndContour = surfc(Axes, X1Plot, X2Plot, reshape(FPred, reverse(Options.NumPlotGrid)), ...
                            'FaceColor', 'r', 'FaceAlpha', .5, 'LineStyle', '-');
                        Curves = [Curves, {SurfAndContour(1), SurfAndContour(2)}];
                        % (6) Plot next point
                        XNextToPlot = makeXTablePlottable(Results, Results.NextPoint);
                        MarkerSpec = {'o', 'MarkerSize', 8, 'Color', 'black', ...
                            'MarkerFaceColor', 'black'};
                        Curves{end+1} = plot3(Axes, XNextToPlot(1), XNextToPlot(2), predictObjective(Results, Results.NextPoint), ...
                            MarkerSpec{:});
                        % Plot model min if exists
                        if ~isempty(ModelMinXTable)
                            XMinPlottable = makeXTablePlottable(Results, ModelMinXTable);
                            Curves{end+1} = plot3(Axes, XMinPlottable(1), XMinPlottable(2), ModelMinF, 'r*');
                        end
                        % (7) Label axes
                        xlabel(Axes, VarSpec.Names{1});
                        ylabel(Axes, VarSpec.Names{2});
                        % (8) set axes scale and limits
                        if isequal(VarSpec.Transforms{1}, 'log')
                            Axes.XScale = 'log';
                            Axes.XLim = [VarSpec.LBs(1), VarSpec.UBs(1)];
                        end
                        if isequal(VarSpec.Transforms{2}, 'log')
                            Axes.YScale = 'log';
                            Axes.YLim = [VarSpec.LBs(2), VarSpec.UBs(2)];
                        end
                        if isequal(VarSpec.Types{1}, 'categorical')
                            cats = categories(VarSpec.Categories{1});
                            Axes.XTick = 1:numel(cats);
                            Axes.XTickLabel = cats;
                            Axes.XTickLabelRotation = 45;
                        end
                        if isequal(VarSpec.Types{2}, 'categorical')
                            cats = categories(VarSpec.Categories{2});
                            Axes.YTick = 1:numel(cats);
                            Axes.YTickLabel = cats;
                            Axes.YTickLabelRotation = 45;
                        end
                        
                end
            end
            
            
        end
    end
        
    %% Internal
    properties(Hidden, SetAccess=protected)
        % Model training data
        XTrain;
        FTrain;
        ObjectiveEvaluationTimeTrain;
        ConstraintTrain;
        ErrorTrain;
        QueuedX;
        % Models
        ObjectiveFcnGP;
        ObjectiveEvaluationTimeGP;
        ConstraintGPs;
        ErrorGP;
        % Current iteration info
        IncumbentF;
        IncumbentX;
        XNext;
        PlotFcnStop;
        OutputFcnStop
        % other
        SampledGridIndices;
        IterationStartTime;
        ObjectiveEvaluationTime;
        PrivOptions;
        PrivIterationTimeTrace;
        PrivIndexOfMinimumTrace;
        PrivFMinTrace;
        PrivFMinEstimatedTrace;
        PrivUserDataTrace;
    end

    methods     % Dependent property getters
        function ObjectiveFcn = get.ObjectiveFcn(this)
            ObjectiveFcn = this.PrivOptions.ObjectiveFcn;
        end
        
        function VariableDescriptions = get.VariableDescriptions(this)
            VariableDescriptions = this.PrivOptions.VariableDescriptions;
        end
        
        function this = set.VariableDescriptions(this, VariableDescriptions)
            this.PrivOptions.VariableDescriptions = VariableDescriptions;
        end
        
        function s = get.Options(this)
            % Return all the visible public options properties in a struct.
            props = properties(this.PrivOptions);
            for i = 1:numel(props)
                s.(props{i}) = this.PrivOptions.(props{i});
            end
        end
        
        function XAtMinObjective = get.XAtMinObjective(this)
            XAtMinObjective = bestPoint(this, 'Criterion', 'min-observed');
        end
        
        function MinObjective = get.MinObjective(this)
            [~, MinObjective] = bestPoint(this, 'Criterion', 'min-observed');
        end
        
        function XAtMinEstimatedObjective = get.XAtMinEstimatedObjective(this)
            XAtMinEstimatedObjective = bestPoint(this);
        end
        
        function MinEstimatedObjective = get.MinEstimatedObjective(this)
            MinEstimatedObjective = predictObjective(this, this.XAtMinEstimatedObjective);
        end
        
        function NumObjectiveEvaluations = get.NumObjectiveEvaluations(this)
            NumObjectiveEvaluations = size(this.XTrain,1);
        end
        
        function TotalElapsedTime = get.TotalElapsedTime(this)
            TotalElapsedTime = nansum(this.IterationTimeTrace);
        end
        
        function NextPoint = get.NextPoint(this)
            if isempty(this.XNext)
                NextPoint = table;
            else
                NextPoint = untransformPoints(this, this.XNext, true);
                NextPoint = applyConditionalVariableFcn(this, NextPoint);
                NextPoint = canonicalizePoints(this, NextPoint);
            end
        end
        
        function XTrace = get.XTrace(this)
            XTrace = conditionalizeX(this, this.XTrain);
        end
        
        function ObjectiveTrace = get.ObjectiveTrace(this)
            ObjectiveTrace = this.FTrain(:);
        end
        
        function ObjectiveEvaluationTimeTrace = get.ObjectiveEvaluationTimeTrace(this)
            ObjectiveEvaluationTimeTrace = this.ObjectiveEvaluationTimeTrain(:);
        end
        
        function ConstraintsTrace = get.ConstraintsTrace(this)
            ConstraintsTrace = this.ConstraintTrain;
        end
        
        function UserDataTrace = get.UserDataTrace(this)
            UserDataTrace = this.PrivUserDataTrace;
        end
        
        function ErrorTrace = get.ErrorTrace(this)
            ErrorTrace = this.ErrorTrain;
        end
        
        function FeasibilityTrace = get.FeasibilityTrace(this)
            FeasibilityTrace = allConstraintsSatisfied(this, this.XTrain, ...
                this.ConstraintGPs, this.ErrorGP);
        end
        
        function FeasibilityProbabilityTrace = get.FeasibilityProbabilityTrace(this)
            FeasibilityProbabilityTrace(:,1) = ProbAllConstraintsSatisfied(this, this.XTrain);
        end
        
        function IndexOfMinimumTrace = get.IndexOfMinimumTrace(this)
            IndexOfMinimumTrace = this.PrivIndexOfMinimumTrace(:);
        end
        
        function ObjectiveMinimumTrace = get.ObjectiveMinimumTrace(this)
            ObjectiveMinimumTrace = this.PrivFMinTrace(:);
        end
        
        function EstimatedObjectiveMinimumTrace = get.EstimatedObjectiveMinimumTrace(this)
            EstimatedObjectiveMinimumTrace = this.PrivFMinEstimatedTrace(:);
        end
        
        function IterationTimeTrace = get.IterationTimeTrace(this)
            IterationTimeTrace = this.PrivIterationTimeTrace(:);
        end
    end
    
    methods(Access=protected)
        %% Main algorithm
        function this = run(this)
            checkForOptimizableVariables(this);
            checkXConstraintFcnSatisfiability(this);
            this.IterationStartTime = tic;
            this = processInitializationData(this);
            this = callPlotFcn(this, 'initial');
            this = callOutputFcn(this, 'initial');
            this = fitModels(this);
            [this.IncumbentF, this.IncumbentX] = findIncumbent(this);
            this = chooseNextPoint(this);
            iteration = this.NumObjectiveEvaluations + 1;
            
            this.customOpts.trajectory.data = cell(0);
            this.customOpts.trajectory.policy = [];
            
            while ~optimizationFinished(this, iteration)
                % Perform feval
                this = performFcnEval(this);
                
                %this.customOpts.trajectory = this.PrivUserDataTrace;
                
                this.customOpts.trajectory.data(end+1,:) = this.PrivUserDataTrace{end,1}.data(1,:);
                this.customOpts.trajectory.policy(end+1,:) = this.PrivUserDataTrace{end,1}.policy(1,:);
                
%                 if this.MinObjective == this.ObjectiveTrace(end)
%                     visCartPole(this.PrivUserDataTrace{end,1}.data{1,1}, this.customOpts.bounds);
%                 end
                
                this = fitModels(this);
                % Update Objective minimum traces
                [~, MinObsObjective, MinObsLoc] = bestPoint(this, 'Criterion', 'min-observed');
                this.PrivIndexOfMinimumTrace(iteration) = MinObsLoc;
                this.PrivFMinTrace(iteration) = MinObsObjective;
                this.PrivFMinEstimatedTrace(iteration) = this.MinEstimatedObjective;
                % Choose next point
                [this.IncumbentF, this.IncumbentX] = findIncumbent(this);
                this = chooseNextPoint(this);
                % Update timing, plots and output
                this.PrivIterationTimeTrace(iteration) = toc(this.IterationStartTime);
                this.IterationStartTime = tic;
                printProgress(this, iteration);
                this = callPlotFcn(this, 'iteration');
                this = callOutputFcn(this, 'iteration');
                iteration = iteration + 1;
            end
            this = pruneGridIndices(this);
            this = callPlotFcn(this, 'done');
            this = callOutputFcn(this, 'done');
            showStoppingReason(this);
            showBestPoints(this);
        end
        
        function this = processInitializationData(this)
            % Set this.QueuedX to be a set of points to evaluate first.
            if this.NumObjectiveEvaluations > 0
                % Only process initialization data if not resuming.
                return;
            end
            Opts = this.PrivOptions;
            if isempty(Opts.InitialX)
                % No initial points passed. Choose initial points.
                if isequal(Opts.AcquisitionFunctionName, 'grid')
                    % Random grid points
                    [X, this] = gridXFeasiblePoint(this);
                    XPoints = X;
                    while ~isempty(X) && size(XPoints,1) < Opts.NumSeedPoints
                        [X, this] = gridXFeasiblePoint(this);
                        if ~isempty(X)
                            XPoints(end+1, :) = X;
                        end
                    end
                else
                    % Random initial points
                    XPoints = initialXFeasiblePoints(this, this.PrivOptions.NumSeedPoints);
                end
                this.QueuedX = XPoints;
            elseif isempty(Opts.InitialObjective)
                % Initial X passed only. Recall that it is a table.
                this.QueuedX = transformPoints(this, Opts.InitialX);
            else
                % Initial evaluated points passed. use them directly.
                XTbl = canonicalizePoints(this, Opts.InitialX);
                this.XTrain = transformPoints(this, XTbl);
                this.FTrain = Opts.InitialObjective;
                this.ObjectiveEvaluationTimeTrain = Opts.InitialObjectiveEvaluationTimes;
                this.PrivIndexOfMinimumTrace(1:this.NumObjectiveEvaluations,1) = NaN;
                this.PrivFMinTrace(1:this.NumObjectiveEvaluations,1) = NaN;
                this.PrivFMinEstimatedTrace(1:this.NumObjectiveEvaluations,1) = NaN;
                this.PrivIterationTimeTrace = Opts.InitialIterationTimes;
                this.PrivUserDataTrace = Opts.InitialUserData;
                this.ConstraintTrain = Opts.InitialConstraintViolations;
                this.ErrorTrain = Opts.InitialErrorValues;
            end
        end
        
        function XBest = initialXFeasiblePoints(this, N)
            % Use random search to find N X-feasible points within bounds
            % for which the smallest Euclidean distance between points is
            % maximized.
            LB = this.PrivOptions.VarSpec.LBTrans;
            UB = this.PrivOptions.VarSpec.UBTrans;
            % Random search
            XBest = randomXFeasiblePoints(this, N);
            bsf = minDist(XBest);
            reps = 50;
            for r = 1:reps
                X = randomXFeasiblePoints(this, N);
                if minDist(X) > bsf
                    bsf = minDist(X);
                    XBest = X;
                end
            end
            function d = minDist(X)
                X = (X-LB)./(UB-LB);
                PointDists = pdist(X);
                d = min(PointDists(:));
            end
        end

        function this = chooseNextPoint(this)
            % Choose the next point to evaluate
            if this.NumObjectiveEvaluations < size(this.QueuedX,1)
                this.XNext = this.QueuedX(this.NumObjectiveEvaluations + 1, :);
            else
                [X, ChoseRandom, this] = chooseNextPointInternal(this);
                % Perform overexploitation loop. Overexploitation is impossible
                % if we have no IncumbentF, because we haven't found any
                % feasible points yet.
                if ~isnan(this.IncumbentF) && ~ChoseRandom && ...
                        ismember(this.PrivOptions.AcquisitionFunctionName,{'expectedimprovementplus','expectedimprovementpersecondplus'}) && ...
                        exploitingTooMuch(this, X)
                    % Save ObjectiveFcnGP
                    ObjFcnGP = this.ObjectiveFcnGP;
                    iteration = 1;
                    BullAdjustmentN = this.NumObjectiveEvaluations;
                    while exploitingTooMuch(this, X) && iteration <= this.PrivOptions.MaxExploitIterations
                        if this.PrivOptions.Verbose >= 2
                            fprintf('\n');
                            bayesoptim.printInfo('Overexploit');
                            disp(conditionalizeX(this, X));
                        end
                        % Refit GPs
                        this = fitGPModels(this, true, BullAdjustmentN);
                        % Choose the next point (maximize Acquisition Function)
                        [X, ~, this] = chooseNextPointInternal(this);
                        % Update loop
                        iteration = iteration + 1;
                        BullAdjustmentN = BullAdjustmentN*10;
                    end
                    % Restore ObjectiveFcnGP
                    this.ObjectiveFcnGP = ObjFcnGP;
                end
                this.XNext = X;
            end
        end
        
        function [tf, ReasonKey] = shouldChooseRandomPoint(this)
            tf = false;
            ReasonKey = [];
            if isequal(this.PrivOptions.AcquisitionFunctionName, 'grid')
                return;
            elseif isequal(this.PrivOptions.AcquisitionFunctionName, 'random')
                % AF is 'random'
                tf = true;
                ReasonKey = 'RandomPointRandAcq';
            elseif nansum(this.ErrorTrain < 0) < this.PrivOptions.NumSeedPoints
                % Not enough non-error points
                tf = true;
                ReasonKey = 'RandomPointInsuffSeedPoints';
            elseif isempty(this.ObjectiveFcnGP)
                % No Objective model
                tf = true;
                ReasonKey = 'RandomPointNoObjModel';
            else
                SigmaF = iGetSigmaF(this.ObjectiveFcnGP);
                if SigmaF/(SigmaF + this.ObjectiveFcnGP.Sigma) < this.PrivOptions.SigmaFTol
                    % Flat Objective model
                    tf = true;
                    ReasonKey = 'RandomPointFlatObjModel';
                end
            end
        end
        
        function [XBest, ChooseRandom, this] = chooseNextPointInternal(this)
            [ChooseRandom, ReasonKey] = shouldChooseRandomPoint(this);
            if ChooseRandom
                if this.PrivOptions.Verbose >= 2
                    bayesoptim.printInfo('ChoosingRandomPoint');
                    bayesoptim.printInfo(ReasonKey);
                end
                XBest = randomXFeasiblePoint(this);
            elseif isequal(this.PrivOptions.AcquisitionFunctionName, 'grid')
                [XBest, this] = gridXFeasiblePoint(this);
            else
                % Choose the next point by maximizing the
                % constraint-weighted acquisition function.
                switch this.PrivOptions.AcquisitionFunctionName
                    case 'probabilityofimprovement'
                        AFcn = @(X)bayesoptim.probabilityOfImprovement(X, this.ObjectiveFcnGP, ...
                            this.IncumbentF, this.ObjectiveFcnGP.Sigma);
                    case {'expectedimprovement', 'expectedimprovementplus'}
                        AFcn = @(X)bayesoptim.expectedImprovement(X, this.ObjectiveFcnGP, ...
                            this.IncumbentF);
                    case 'lowerconfidencebound'
                        Kappa = 2;
                        AFcn = @(X)bayesoptim.lowerConfidenceBound(X, this.ObjectiveFcnGP, Kappa);
                    case {'expectedimprovementpersecond', 'expectedimprovementpersecondplus'}
                        AFcn = @(X)bayesoptim.expectedImprovementPerSecond(X, this.ObjectiveFcnGP, ...
                            this.ObjectiveEvaluationTimeGP, this.IncumbentF);
                    otherwise
                        bayesoptim.err('AFUnknown', this.PrivOptions.AcquisitionFunctionName);
                end
                VarSpec = this.PrivOptions.VarSpec;
                XBest = iFminbndGlobal(@constraintWeightedNegAF, VarSpec.LBTrans, ...
                    VarSpec.UBTrans, this.PrivOptions.NumRestartCandidates, ...
                    this.PrivOptions.NumRestarts, this.PrivOptions.VerboseRestarts, ...
                    this.PrivOptions.MaxIterPerRestart, this.PrivOptions.RelTol);
            end
            if ~isempty(XBest)
                XBest = legalizePoints(this, XBest);
            end
            
            function NegAF = constraintWeightedNegAF(X)
                % Return the constraint-weighted negative AF value for all
                % feasible points in X (to be minimized). NegAF(i) = Inf if
                % not feasible
                NegAF = inf(size(X,1),1);
                [Xcanon, InputFeasible] = legalizePoints(this, X);
                if any(InputFeasible)
                    Xfeasible = Xcanon(InputFeasible,:);
                    if isnan(this.IncumbentF)
                        % We have constraints, but no incumbent, meaning that
                        % no feasible points were found. Maximize only the
                        % constraint probability.
                        NegAF(InputFeasible) = -ProbAllConstraintsSatisfied(this, Xfeasible);
                    else
                        % Constraints and an incumbent. Maximize constraint-weighted AF.
                        NegAF(InputFeasible) = -AFcn(Xfeasible) .*  ProbAllConstraintsSatisfied(this, Xfeasible);
                    end
                end
            end
        end
        
        function [XBest, this] = sampleGridWithoutReplacement(this)
            % Return [] if the grid has already been fully sampled
            Divs = this.PrivOptions.NumGridDivisions;
            if size(this.SampledGridIndices, 1) == prod(Divs)
                if this.PrivOptions.Verbose >= 2
                    bayesoptim.printInfo('GridFullySampled');
                end
                this.SampledGridIndices = [];
            end
            Indices = iFindUnusedGridIndices(Divs, this.SampledGridIndices);
            XBest = XFromGridIndices(this, Indices);
            this.SampledGridIndices(end+1, :) = Indices;
        end
        
        function this = pruneGridIndices(this)
            if isequal(this.PrivOptions.AcquisitionFunctionName, 'grid') && ~isempty(this.XNext)
                this.SampledGridIndices(end, :) = [];
            end
        end
        
        function X = XFromGridIndices(this, Indices)
            % Grid points include both endpoints. Indices are origin 1.
            GridRes = this.PrivOptions.NumGridDivisions;
            LB = this.PrivOptions.VarSpec.LBTrans;
            UB = this.PrivOptions.VarSpec.UBTrans;
            X = (Indices-1)./(GridRes-1) .* (UB-LB) + LB;
        end
        
        function TF = exploitingTooMuch(this, X)
            % X is a design matrix. TF is a logical vector. TF(i)=true if
            % the function std is less than alpha times the noise std at
            % the point X(i,:).
            %
            % "You can't overexploit an infeasible point": If X is being
            % evaluated for overexploitation, then it has been chosen by
            % the AF as the most desirable point. If this point is not
            % known to be feasible, then there is something to be gained by
            % evaluating it. We need to know whether it is feasible or not.
            if isempty(this.ObjectiveFcnGP)
                TF = false(size(X,1),1);
            else
                % Check ObjectiveFcnGP condition
                [~, YSD] = predict(this.ObjectiveFcnGP, X);
                FSD = funStd(YSD, this.ObjectiveFcnGP.Sigma);
                TF = FSD < this.PrivOptions.ExplorationRatio*this.ObjectiveFcnGP.Sigma;
                
                % Check feasibility:
                TF = TF && allConstraintsSatisfied(this, X, this.ConstraintGPs, this.ErrorGP);
            end
        end
        
        function tf = optimizationFinished(this, iteration)
            if iteration == 1
                tf = false;
            elseif isempty(this.XNext)
                tf = true;
            elseif this.PlotFcnStop || this.OutputFcnStop || ...
                    this.NumObjectiveEvaluations >= this.PrivOptions.MaxObjectiveEvaluations || ...
                    this.TotalElapsedTime >= this.PrivOptions.MaxTime
                tf = true;
            else
                tf = false;
            end
        end
        
        %% Function evaluation
        function this = performFcnEval(this)
            % Perform a function evaluation on XNext and store the results
            Opts = this.PrivOptions;
            if Opts.Verbose >= 2
                fprintf('__________________________________________________________\n');
                bayesoptim.printInfo('IterationMessage', length(this.FTrain)+1);
                disp(conditionalizeX(this, this.XNext));
            end
            % Call fcn
            [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ...
                ObjectiveFcnObjectiveEvaluationTime, this] = callObjFcn(this, this.XNext);
            % Store results in traces
            this.XTrain(end+1, :) = legalizePoints(this, this.XNext);
            this.FTrain(end+1) = Objective;
            if this.PrivOptions.NumCoupledConstraints > 0
                this.ConstraintTrain(end+1,:) = ConstraintViolations;
            end
            this.PrivUserDataTrace{end+1,1} = UserData;
            this.ErrorTrain(end+1,1) = ErrorConstraintViolation;
            this.ObjectiveEvaluationTimeTrain(end+1) = ObjectiveFcnObjectiveEvaluationTime;
        end
        
        function [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ...
                ObjectiveEvaluationTime, this] = callObjFcn(this, X)
            if ~isempty(this.PrivOptions.ObjectiveNargout)
                % nargout is known. Call ObjectiveFcn normally.
                [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ObjectiveEvaluationTime, this] ...
                    = callObjNormally(this, X);
            else
                % nargout is unknown. Try nargout=3:
                StartTic = tic;
                ErrorConstraintViolation = -1;
                try
                    [Objective, ConstraintViolations, UserData] = this.ObjectiveFcn(conditionalizeX(this, X));
                catch msg
                    if ismember(msg.identifier, {'MATLAB:maxlhs', 'MATLAB:unassignedOutputs', 'MATLAB:TooManyOutputs'})
                        % nargout < 3.
                        % Warn the user if we had to waste a function evaluation
                        if ismember(msg.identifier, {'MATLAB:maxlhs', 'MATLAB:unassignedOutputs'})
                            bayesoptim.warn('UnknownNargoutWasExpensive');
                        end
                        % Set nargout based on NumCoupledConstraints, call recursively, and return:
                        this.PrivOptions.ObjectiveNargout = nargoutFromNumConstraints(this.PrivOptions.NumCoupledConstraints);
                        [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ObjectiveEvaluationTime, this] ...
                            = callObjFcn(this, X);
                        return;
                    else
                        rethrow(msg);
                    end
                end
                % Calling with nargout=3 succeeded. Set nargout and finish up.
                this.PrivOptions.ObjectiveNargout = 3;
                [Objective, ConstraintViolations, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
                    = finishObjEval(this, Objective, ConstraintViolations, -1, StartTic);
            end
        end
        
        function [Objective, ConstraintViolations, UserData, ErrorConstraintViolation, ObjectiveEvaluationTime, this] ...
                = callObjNormally(this, X)
            % Call the objective fcn with the correct nargout, catching
            % errors and timing runtime
            StartTic = tic;
            ErrorConstraintViolation = -1;
            try
                switch this.PrivOptions.ObjectiveNargout
                    case 1
                        Objective = this.ObjectiveFcn(conditionalizeX(this, X));
                        ConstraintViolations = [];
                        UserData = [];
                    case 2
                        [Objective, ConstraintViolations] = this.ObjectiveFcn(conditionalizeX(this, X));
                        UserData = [];
                    otherwise
                        [Objective, ConstraintViolations, UserData] = this.ObjectiveFcn(conditionalizeX(this, X));
                end
            catch msg
                if isequal(msg.identifier, 'MATLAB:unassignedOutputs') && this.PrivOptions.NumCoupledConstraints > 0
                    bayesoptim.err('ObjectiveNargoutWrong', this.PrivOptions.NumCoupledConstraints);
                else
                    rethrow(msg);
                end
            end
            [Objective, ConstraintViolations, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
                = finishObjEval(this, Objective, ConstraintViolations, ErrorConstraintViolation, StartTic);
        end
        
        function [Objective, ConstraintViolations, ErrorConstraintViolation, ObjectiveEvaluationTime] ...
                = finishObjEval(this, Objective, ConstraintViolations, ErrorConstraintViolation, StartTic)
            % Check constraint violation dimensions
            if numel(ConstraintViolations) ~= this.PrivOptions.NumCoupledConstraints
                bayesoptim.err('ObjectiveNumConstraintsWrong',...
                    numel(ConstraintViolations), this.PrivOptions.NumCoupledConstraints);
            end
            % Set illegal Objective and ConstraintViolation to NaN
            Objective = iNanIfBad(Objective);
            ConstraintViolations = arrayfun(@iNanIfBad, ConstraintViolations);
            if isnan(Objective)
                ErrorConstraintViolation = 1;
            end
            % Record runtime
            ObjectiveEvaluationTime = toc(StartTic);
        end
        
        %% Modeling
        function this = fitModels(this)
            this = fitGPModels(this, false, 1);
        end
        
        function this = fitGPModels(this, BullAdjustment, BullAdjustmentN)
            if this.PrivOptions.FitModels
                this.ObjectiveFcnGP = fitObjectiveFcnGP(this, BullAdjustment, BullAdjustmentN);
                this.ObjectiveEvaluationTimeGP = fitObjectiveEvaluationTimeGP(this);
                this.ConstraintGPs = fitConstraintGPs(this);
                this.ErrorGP = fitErrorGP(this);
            else
                this.ObjectiveFcnGP = [];
                this.ObjectiveEvaluationTimeGP = [];
                this.ConstraintGPs = cell(1, this.PrivOptions.NumCoupledConstraints);
                this.ErrorGP = [];
            end
        end
        
        function ObjectiveFcnGP = fitObjectiveFcnGP(this, BullAdjustment, BullAdjustmentN)
            if all(isnan(this.FTrain))
                ObjectiveFcnGP = [];
                if this.PrivOptions.Verbose >= 2
                    bayesoptim.printInfo('CantFitObjective');
                end
            else
                [KernelParams0, ConstantKernelParameters] = kernelParamsForFitrgp(this, this.FTrain);
                if this.PrivOptions.IsObjectiveDeterministic
                    SigmaLowerBound = max(1e-8, 1e-4*nanstd(this.FTrain));
                    sigma0 = SigmaLowerBound;
                    ConstantFSigma = true;
                    ObjectiveFcnGP = iFitrgpRobust(this.customOpts,this.XTrain, this.FTrain, SigmaLowerBound, ...
                        'BasisFunction', 'constant', ...
                        'ConstantKernelParameters', ConstantKernelParameters, ...
                        'ConstantSigma', ConstantFSigma,...
                        'FitMethod', 'exact', ...
                        'KernelFunction', 'ardmatern52', ...
                        'KernelParameters', KernelParams0, ...
                        'PredictMethod', 'exact', ...
                        'Sigma', sigma0, ...
                        'Standardize', false,...
                        'Verbose', double(this.PrivOptions.Verbose >= 3));
                else
                    SigmaLowerBound = max(1e-8, nanstd(this.FTrain)*1e-2);
                    sigma0 = max(SigmaLowerBound, nanstd(this.FTrain)/this.PrivOptions.Sigma0Divisor);
                    ObjectiveFcnGP = iFitrgpRobust(this.customOpts,this.XTrain, this.FTrain, SigmaLowerBound, ...
                        'BasisFunction', 'constant', ...
                        'ConstantKernelParameters', ConstantKernelParameters, ...
                        'FitMethod', 'exact', ...
                        'KernelFunction', 'ardmatern52', ...
                        'KernelParameters', KernelParams0, ...
                        'PredictMethod', 'exact', ...
                        'Sigma', sigma0, ...
                        'Standardize', false,...
                        'Verbose', double(this.PrivOptions.Verbose >= 3));
                end
                if BullAdjustment
                    % Apply Bull's adjustment to the global kernel amplitude, and
                    % recompute alphas. This requries passing in all the current values
                    % of other parameters and using 'FitMethod' 'none'.
                    KernelParameters =  ObjectiveFcnGP.KernelInformation.KernelParameters;
                    KernelParameters(end) = KernelParameters(end)*BullAdjustmentN;
                    % recompute alphas only
                    ObjectiveFcnGP = iFitrgpRobust(this.customOpts,this.XTrain, this.FTrain, SigmaLowerBound,...
                        'KernelFunction', 'ardmatern52', ...
                        'KernelParameters', KernelParameters, ...
                        'Sigma', ObjectiveFcnGP.Sigma, ...
                        'BasisFunction', 'constant', ...
                        'Beta', ObjectiveFcnGP.Beta, ...
                        'FitMethod', 'none', ...
                        'PredictMethod', 'exact', ...
                        'Standardize', false);
                end
            end
        end
        
        function ObjectiveEvaluationTimeGP = fitObjectiveEvaluationTimeGP(this)
            % Only fit model to successful function evaluations
            FitIndices = ~isnan(this.ObjectiveEvaluationTimeTrain) & ~isnan(this.FTrain);
            if ~any(FitIndices)
                ObjectiveEvaluationTimeGP = [];
                if this.PrivOptions.Verbose >= 2
                    bayesoptim.printInfo('CantFitObjectiveEvaluationTime');
                end
            else
                [KernelParams0, ConstantKernelParameters] = kernelParamsForFitrgp(this, this.ObjectiveEvaluationTimeTrain);
                SigmaLowerBound = std(log(this.ObjectiveEvaluationTimeTrain(FitIndices)))*1e-2;
                ObjectiveEvaluationTimeGP = iFitrgpRobust(this.customOpts,this.XTrain(FitIndices,:), ...
                    log(this.ObjectiveEvaluationTimeTrain(FitIndices)), ...
                    SigmaLowerBound, ...
                    'BasisFunction', 'constant', ...
                    'ConstantKernelParameters', ConstantKernelParameters, ...
                    'FitMethod', 'exact', ...
                    'KernelFunction', 'ardmatern52', ...
                    'KernelParameters', KernelParams0, ...
                    'PredictMethod', 'exact', ...
                    'Standardize', false,...
                    'Verbose', double(this.PrivOptions.Verbose >= 3));
            end
        end
        
        function ConstraintGPs = fitConstraintGPs(this)
            ConstraintGPs = {};
            if ~isempty(this.ConstraintTrain)
                for i = this.PrivOptions.NumCoupledConstraints:-1:1
                    ConstraintGPs{i} = fitConstraintGP(this, this.ConstraintTrain(:,i), ...
                        this.PrivOptions.AreCoupledConstraintsDeterministic(i));
                end
            end
        end
        
        function ConstraintGP = fitConstraintGP(this, ConstraintTrain, isDeterministic)
            if all(isnan(ConstraintTrain))
                ConstraintGP = [];
                if this.PrivOptions.Verbose >= 2
                    bayesoptim.printInfo('CantFitConstraint');
                end
            else
                [KernelParams0, ConstantKernelParameters] = kernelParamsForFitrgp(this, ConstraintTrain);
                if isDeterministic
                    std = nanstd(ConstraintTrain);
                    if std == 0
                        SigmaLowerBound = 1e-4;
                    else
                        SigmaLowerBound = 1e-4*std;
                    end
                    sigma0 = SigmaLowerBound;
                    ConstantSigma = true;
                    ConstraintGP = iFitrgpRobust(this.customOpts,this.XTrain, ConstraintTrain, SigmaLowerBound, ...
                        'BasisFunction', 'constant', ...
                        'ConstantKernelParameters', ConstantKernelParameters, ...
                        'ConstantSigma', ConstantSigma,...
                        'FitMethod', 'exact', ...
                        'KernelFunction', 'ardmatern52', ...
                        'KernelParameters', KernelParams0, ...
                        'PredictMethod', 'exact', ...
                        'Sigma', sigma0, ...
                        'Standardize', false,...
                        'Verbose', double(this.PrivOptions.Verbose >= 3));
                else
                    std = nanstd(ConstraintTrain);
                    if std == 0
                        SigmaLowerBound = 1e-2;
                    else
                        SigmaLowerBound = 1e-2*std;
                    end
                    ConstraintGP = iFitrgpRobust(this.customOpts,this.XTrain, ConstraintTrain, SigmaLowerBound, ...
                        'BasisFunction', 'constant', ...
                        'ConstantKernelParameters', ConstantKernelParameters, ...
                        'FitMethod', 'exact', ...
                        'KernelFunction', 'ardmatern52', ...
                        'KernelParameters', KernelParams0, ...
                        'PredictMethod', 'exact', ...
                        'Standardize', false,...
                        'Verbose', double(this.PrivOptions.Verbose >= 3));
                end
            end
        end
        
        function ErrorGP = fitErrorGP(this)
            if all(isnan(this.ErrorTrain))
                ErrorGP = [];
                if this.PrivOptions.Verbose >= 2
                    bayesoptim.printInfo('CantFitError');
                end
            elseif all(this.ErrorTrain < 0)
                ErrorGP = [];
            else
                [KernelParams0, ConstantKernelParameters] = kernelParamsForFitrgp(this, this.ErrorTrain);
                SigmaLowerBound = nanstd(this.ErrorTrain)*1e-4;
                ErrorGP = iFitrgpRobust(this.customOpts,this.XTrain, this.ErrorTrain, SigmaLowerBound, ...
                    'BasisFunction', 'constant', ...
                    'ConstantKernelParameters', ConstantKernelParameters, ...
                    'FitMethod', 'exact', ...
                    'KernelFunction', 'ardmatern52', ...
                    'KernelParameters', KernelParams0, ...
                    'PredictMethod', 'exact', ...
                    'Standardize', false,...
                    'Verbose', double(this.PrivOptions.Verbose >= 3));
            end
        end
        
        %% Proxy optimization
        function [IncumbentF, IncumbentX] = findIncumbent(this, BOOptions)
            % Find the point that minimizes the mean of the Objective
            % function model. Optionally accept modified Options.
            if nargin < 2
                BOOptions = this.PrivOptions;
            end
            if isempty(this.ObjectiveFcnGP) || ...
                    (this.PrivOptions.NumCoupledConstraints > 0 && any(cellfun(@isempty, this.ConstraintGPs)))
                IncumbentF = NaN;
                IncumbentX = [];
            else
                VarSpec = this.PrivOptions.VarSpec;
                IncumbentX = iFminbndGlobal(@(X)GPFMeanOnPoints(this,X), VarSpec.LBTrans, VarSpec.UBTrans, ...
                    BOOptions.NumRestartCandidates, BOOptions.NumRestarts, BOOptions.VerboseRestarts, ...
                    BOOptions.MaxIterPerRestart, BOOptions.RelTol);
                % XAtMin is not discretized (although it was discretized during optimization)
                IncumbentX = legalizePoints(this, IncumbentX);
                IncumbentF = GPFMeanOnPoints(this, IncumbentX);
                % If no finite incumbent was found, set it to NaN
                if ~isfinite(IncumbentF)
                    IncumbentF = NaN;
                    IncumbentX = [];
                end
            end
        end
        
        function Objective = GPFMeanOnPoints(this, X)
            % Return the mean ObjectiveFcnGP value for all feasible points
            % in X. Return Inf for infeasible points. Note that we must
            % transform to native space in order to apply the user-supplied
            % conditionalVariableFcn and constraint functions, and after
            % that must canonicalize X and transform back before we can
            % apply the GP models.
            [Xcanon, InputFeasible] = legalizePoints(this, X);
            allSat = allConstraintsSatisfied(this, Xcanon, this.ConstraintGPs, this.ErrorGP);
            feasible = InputFeasible & allSat;
            Objective = inf(size(X,1),1);
            Objective(feasible) = predict(this.ObjectiveFcnGP, Xcanon(feasible,:));
        end

        function TF = isInputFeasible(this, X)
            % TF(i) is true iff X(i,:) satisfies the bounds and XConstraintFcn.
            TF = all(bsxfun(@ge, X, this.PrivOptions.VarSpec.LBTrans), 2) & ...
                all(bsxfun(@le, X, this.PrivOptions.VarSpec.UBTrans), 2);
            TF = TF & applyXConstraint(this, untransformPoints(this, X, true));
        end
        
        function InputFeasible = applyXConstraint(this, XTable)
            if isempty(this.PrivOptions.XConstraintFcn)
                InputFeasible = true(height(XTable),1);
            else
                try
                    InputFeasible = this.PrivOptions.XConstraintFcn(XTable);
                catch me
                    bayesoptim.warn('XConstraintFcnError');
                    rethrow(me);
                end
                if ~(islogical(InputFeasible) && isequal(size(InputFeasible), [height(XTable),1]))
                    bayesoptim.err('XConstraintFcnOutput', num2str(height(XTable)), class(InputFeasible), num2str(size(InputFeasible)));
                end
            end
        end
        
        function checkForOptimizableVariables(this)
            if this.PrivOptions.VarSpec.NumVars == 0
                bayesoptim.err('NoOptimizableVariables');
            end
        end
        
        function checkXConstraintFcnSatisfiability(this)
            % Check satisfiability of the XConstraint
            N = 10000;
            X = iUniformPointsWithinBounds(N, this.PrivOptions.VarSpec.LBTrans, ...
                this.PrivOptions.VarSpec.UBTrans);
            NumFeas = sum(isInputFeasible(this, X));
            if NumFeas == 0
                bayesoptim.err('NoFeasiblePoints', N);
            elseif NumFeas < .01*N
                bayesoptim.warn('FewFeasiblePoints', N, NumFeas);
            end
        end
        
        function TF = allConstraintsSatisfied(this, Xcanon, ConstraintGPs, ErrorGP)
            % Return true for every row of Xcanon that satisfies all of the
            % constraints to their required tolerances. Empty constraints
            % are never satisfied. XCanon is canonicalized: It is
            % well-formed and has NaNs filled with LBs.
            if numel(ConstraintGPs) == 0
                Sat = true;
            else
                for i = numel(ConstraintGPs):-1:1
                    if isempty(ConstraintGPs{i})
                        Sat(1:size(Xcanon,1), i) = false;
                    else
                        [means, ysds] = predict(ConstraintGPs{i}, Xcanon);
                        fsds = funStd(ysds, ConstraintGPs{i}.Sigma);
                        Sat(:,i) = normcdf(0, means, fsds) >= 1-this.PrivOptions.CoupledConstraintTolerances(i);
                    end
                end
            end
            if isempty(ErrorGP)
                % Assume no error if no error GP
                ErrorConstraintSat(1:size(Xcanon,1), 1) = true;
            else
                [means, ysds] = predict(ErrorGP, Xcanon);
                fsds = funStd(ysds, ErrorGP.Sigma);
                ErrorConstraintSat = normcdf(0, means, fsds) >= 1-this.PrivOptions.ErrorConstraintTol;
            end
            TF = all(Sat, 2) & ErrorConstraintSat;
        end
        
        %% Output and plotting
        % Plotting helper methods
        function printProgress(this, iteration)
            if this.PrivOptions.Verbose >= 1
                % Compute Action
                [~,~,MinObsLoc] = bestPoint(this, 'Criterion', 'min-observed');
                if isnan(this.FTrain(end))
                    Action = 'Error ';
                elseif ~isempty(this.ConstraintsTrace) && any(this.ConstraintsTrace(end,:) > 0)
                    Action = 'Infeas';
                elseif MinObsLoc == iteration
                    Action = 'Best  ';
                else
                    Action = 'Accept';
                end
                XTable = this.XTrace(end,:);
                
                % Display header.
                IterEvalWidth = 7+9;
                VarNameFieldWidth = 12;
                IncludeEstimCol = ~ismember(this.PrivOptions.AcquisitionFunctionName, {'grid','random'});
                MandatoryFields = 3+IncludeEstimCol;
                if rem(iteration, this.PrivOptions.DisplayHeaderInterval) == 1
                    % Header top bar
                    NumEquals = IterEvalWidth + MandatoryFields*(VarNameFieldWidth+1) - 1 ...
                                + (this.PrivOptions.NumCoupledConstraints + this.NumVars)*(VarNameFieldWidth+3);
                    Bar = [sprintf('|'), sprintf(repmat('=', 1, NumEquals)), sprintf('|\n')];
                    fprintf(Bar);
                    % Header line 1
                    fprintf('| Iter | Eval   | Objective  | Objective  | BestSoFar  |');
                    if IncludeEstimCol
                        fprintf(' BestSoFar  |');
                    end
                    for c = 1:this.PrivOptions.NumCoupledConstraints
                        fprintf(' Constraint%d  |', c);
                    end
                    % Variable names in header line 1
                    for v = 1:this.NumVars
                        Name = this.PrivOptions.VarSpec.Names{v};
                        if numel(Name) > VarNameFieldWidth
                            ControlString = sprintf(' %%%ds-|', VarNameFieldWidth);
                        else
                            ControlString = sprintf(' %%%ds |', VarNameFieldWidth);
                        end
                        fprintf(ControlString, Name(1:min(numel(Name), VarNameFieldWidth)));
                    end
                    fprintf('\n');
                    % Header line 2
                    fprintf('|      | result |            | runtime    | (observed) |');
                    if IncludeEstimCol
                        fprintf(' (estim.)   |');
                    end
                    for c = 1:this.PrivOptions.NumCoupledConstraints
                        fprintf('              |');
                    end
                    % Variable names in header line 2
                    for v = 1:this.NumVars
                        Name = this.PrivOptions.VarSpec.Names{v};
                        if numel(Name) > VarNameFieldWidth
                            str = Name(VarNameFieldWidth+1 : min(numel(Name), 2*VarNameFieldWidth));
                        else
                            str = '';
                        end
                        ControlString = sprintf(' %%-%ds |', VarNameFieldWidth);
                        fprintf(ControlString, str);
                    end
                    fprintf('\n');
                    % Header bottom bar
                    fprintf(Bar);
                end
                
                % Display a line of data
                fprintf('| %4d | %s | %10.5g | %10.5g | %10.5g |', ...
                    iteration, Action, this.FTrain(end), this.ObjectiveEvaluationTimeTrain(end), ...
                    this.MinObjective);
                if IncludeEstimCol
                    fprintf(' %10.5g |', this.MinEstimatedObjective);
                end
                for c = 1:this.PrivOptions.NumCoupledConstraints
                    fprintf(' %12.3g |', this.ConstraintTrain(end,c));
                end
                for v = 1:this.NumVars
                    if iscategorical(XTable.(v))
                        if ismissing(XTable(1,v))
                            Name = '-';
                        else
                            Name = char(XTable.(v));
                        end
                        fprintf(' %12s |', Name(1:min(numel(Name),12)));
                    else
                        if isnan(XTable.(v))
                            Name = '-';
                            fprintf(' %12s |', Name(1:min(numel(Name),12)));
                        else
                            fprintf(' %12.5g |', XTable.(v));
                        end
                    end
                end
                fprintf('\n');
            end
        end
        
        function this = callPlotFcn(this, State)
            PlotFcn = this.PrivOptions.PlotFcn;
            this.PlotFcnStop = false;
            if ~isempty(PlotFcn)
                if ~iscell(PlotFcn)
                    PlotFcn = {PlotFcn};
                end
                for i = 1:numel(PlotFcn)
                    try
                        stop = PlotFcn{i}(this, State);
                    catch me
                        bayesoptim.warn('PlotFcnError', i);
                        bayesoptim.printInfo('ErrorDisplay');
                        disp(me.message);
                        stop = false;
                    end
                    if ~(isscalar(stop) && islogical(stop))
                        bayesoptim.warn('PlotFcnOutput', i);
                    end
                    this.PlotFcnStop = this.PlotFcnStop || stop;
                end
                drawnow;
            end
        end
        
        function this = callOutputFcn(this, State)
            OutputFcn = this.PrivOptions.OutputFcn;
            this.OutputFcnStop = false;
            if ~isempty(OutputFcn)
                if ~iscell(OutputFcn)
                    OutputFcn = {OutputFcn};
                end
                for i = 1:numel(OutputFcn)
                    try
                        stop = OutputFcn{i}(this, State);
                    catch me
                        bayesoptim.warn('OutputFcnError', i);
                        bayesoptim.printInfo('ErrorDisplay');
                        disp(me.message);
                        stop = false;
                    end
                    if ~(isscalar(stop) && islogical(stop))
                        bayesoptim.warn('OutputFcnOutput', i);
                    end
                    this.OutputFcnStop = this.OutputFcnStop || stop;
                end
            end
        end
        
        function showStoppingReason(this)
            if this.PrivOptions.Verbose > 0
                fprintf('\n__________________________________________________________\n');
                bayesoptim.printInfo('OptimizationCompleted');
                if this.OutputFcnStop
                    bayesoptim.printInfo('OutputFcnStop');
                end
                if this.PlotFcnStop
                    bayesoptim.printInfo('PlotFcnStop');
                end
                if this.NumObjectiveEvaluations >= this.PrivOptions.MaxObjectiveEvaluations
                    bayesoptim.printInfo('MaxEvalsReached', this.PrivOptions.MaxObjectiveEvaluations);
                end
                if this.TotalElapsedTime >= this.PrivOptions.MaxTime
                    bayesoptim.printInfo('MaxTimeReached', num2str(this.PrivOptions.MaxTime));
                end
                if isempty(this.XNext) && isequal(this.PrivOptions.AcquisitionFunctionName, 'grid')
                    bayesoptim.printInfo('GridSearched');
                end
                bayesoptim.printInfo('TotalEvaluations', this.NumObjectiveEvaluations);
                bayesoptim.printInfo('TotalTime', num2str(this.TotalElapsedTime));
                bayesoptim.printInfo('TotalEvalTime', num2str(nansum(this.ObjectiveEvaluationTimeTrain)));
            end
        end

        function showBestPoints(this)
            if this.PrivOptions.Verbose > 0
                [MinObsXTable, MinObsObjective, MinObsLoc] = bestPoint(this, 'Criterion', 'min-observed');
                if ~isfinite(MinObsLoc)
                    bayesoptim.printInfo('NoFeasibleResult');
                else
                    % Show min observed feasible point
                    fprintf('\n');
                    bayesoptim.printInfo('BestObserved');
                    disp(MinObsXTable);
                    bayesoptim.printInfo('ObservedObjectiveValue', num2str(MinObsObjective));
                    if ~ismember(this.PrivOptions.AcquisitionFunctionName, {'grid','random'})
                        BestPoint = bestPoint(this);
                        if ~isempty(BestPoint)
                            bayesoptim.printInfo('EstimatedObjectiveValue', num2str(predictObjective(this, bestPoint(this))));
                        end
                    end
                    bayesoptim.printInfo('ObservedEvaluationTime', num2str(this.ObjectiveEvaluationTimeTrace(MinObsLoc)));
                    if ~isempty(this.ConstraintsTrace)
                        bayesoptim.printInfoNoReturn('ObservedConstraintViolations');
                        fprintf('[ ');
                        fprintf('%f ', this.ConstraintsTrace(MinObsLoc,:));
                        fprintf(']\n');
                    end
                    % Unless using grid or random, show best point
                    % according to the default criterion
                    if ~ismember(this.PrivOptions.AcquisitionFunctionName, {'grid','random'})
                        BestPoint = bestPoint(this);
                        if ~isempty(BestPoint)
                            fprintf('\n');
                            bayesoptim.printInfo('BestEstimated');
                            disp(BestPoint);
                            bayesoptim.printInfo('EstimatedObjectiveValue', num2str(predictObjective(this, BestPoint)));
                            bayesoptim.printInfo('EstimatedEvaluationTime', num2str(predictObjectiveEvaluationTime(this, BestPoint)));
                            if this.PrivOptions.NumCoupledConstraints > 0
                                ConstraintPred = cellfun(@(GP)predict(GP, transformPoints(this, BestPoint)), ...
                                    this.ConstraintGPs);
                                bayesoptim.printInfoNoReturn('EstimatedConstraintViolations');
                                fprintf('[ ');
                                fprintf('%f ', ConstraintPred);
                                fprintf(']\n');
                            end
                        end
                    end
                    fprintf('\n');
                end
            end
        end
        
        function N = NumVars(this)
            N = this.PrivOptions.VarSpec.NumVars;
        end
        
        function [Xcanon, InputFeasible] = legalizePoints(this, X)
            % Given points X in real, transformed space, return the corresponding
            % conditional-mapped, canonicalized points in the same space.
            XTable = untransformPoints(this, X, true);
            XTable = applyConditionalVariableFcn(this, XTable);
            InputFeasible = applyXConstraint(this, XTable);
            XTable = canonicalizePoints(this, XTable);
            Xcanon = transformPoints(this, XTable);
        end
                
        function XTrain = XTraceToXTrain(this, XTable)
            % Called in resume()
            XTable = applyConditionalVariableFcn(this, XTable);
            XTable = canonicalizePoints(this, XTable);
            XTrain = transformPoints(this, XTable);
        end
        
        function X = makeXTablePlottable(this, XTable)
            % Returns a matrix X that is a plottable version of XTable.
            % NaNs are canonicalized, Categoricals are int-coded.
            XTable = canonicalizePoints(this, XTable);
            for v = 1:width(XTable)
                if iscategorical(XTable.(v))
                    XTable.(v) = intCodeCategorical(XTable.(v));
                end
            end
            X = table2array(XTable);
        end
        
        function XTable = untransformPoints(this, X, HonorInts)
            % Convert a matrix of transformed points to a table in native
            % space. HonorInts=true means that the transformation will
            % round reals to ints (and then to cats if necessary). false
            % means the points are left as reals after transforming to
            % native space. 'false' is only used for plotting meshes.
            if isempty(X)
                XTable = table;
            else
                VarSpec = this.PrivOptions.VarSpec;
                for v = this.NumVars:-1:1
                    switch VarSpec.Types{v}
                        case 'real'
                            if isequal(VarSpec.Transforms{v}, 'log')
                                Vars{v} = exp(X(:,v));
                            else
                                Vars{v} = X(:,v);
                            end
                        case 'integer'
                            if isequal(VarSpec.Transforms{v}, 'log')
                                Vars{v} = exp(X(:,v));
                            else
                                Vars{v} = X(:,v);
                            end
                            if HonorInts
                                % Convert real to int
                                Vars{v} = round(Vars{v});
                            end
                        case 'categorical'
                            Vars{v} = X(:,v);
                            if HonorInts
                                % Convert real to int to cat
                                idx = round(Vars{v});
                                Vars{v} = VarSpec.Categories{v}(idx);
                            end
                    end
                end
                XTable = table(Vars{:}, 'VariableNames', VarSpec.Names);
            end
        end
        
        function Prob = ProbAllConstraintsSatisfied(this, X)
            % The probability that ALL the constraints are satisfied at
            % each row of X. Constraint without models are never satisfied.
            N = size(X,1);
            % Coupled constraints
            if this.PrivOptions.NumCoupledConstraints == 0
                % There are no constraints at all.
                ProbSat = 1;
            else
                for i = numel(this.ConstraintGPs):-1:1
                    if isempty(this.ConstraintGPs{i})
                        % There is a constraint but it has no model yet.
                        ProbSat(1:N,i) = 0;
                    else
                        [means, ysds] = predict(this.ConstraintGPs{i}, X);
                        fsds = funStd(ysds, this.ConstraintGPs{i}.Sigma);
                        ProbSat(:,i) = normcdf(0, means, fsds);
                    end
                end
            end
            % Error constraint
            if isempty(this.ErrorGP)
                % There is no model of error yet.
                ProbErrorSat(1:N,1) = 1;
            else
                [means, ysds] = predict(this.ErrorGP, X);
                fsds = funStd(ysds, this.ErrorGP.Sigma);
                ProbErrorSat(:,1) = normcdf(0, means, fsds);
            end
            % XConstraint
            IsXConstraintSat = satisfiesXConstraint(this, X);
            
            % Calculate prob
            Prob = prod(ProbSat, 2) .* ProbErrorSat .* IsXConstraintSat;
        end
        
        function TF = satisfiesXConstraint(this, X)
            % Returns true for each row for which xConstraint is
            % satisfied.
            XTable = untransformPoints(this, X, true);
            XTable = applyConditionalVariableFcn(this, XTable);
            TF = applyXConstraint(this, XTable);
        end
        
        %% Choosing the best point (used by the bestPoint method)
        function [BestXTable, MinObserved, Iteration] = minObservedPoint(this)
            % Returns the feasible point BestXTable with the minimum
            % observed Objective value. Feasibility is judged by the
            % current constraint and error models.
            FeasMask = double(this.FeasibilityTrace);
            FeasMask(FeasMask==0) = NaN;
            [MinObserved, Iteration] = nanmin(this.ObjectiveTrace .* FeasMask);
            if isfinite(MinObserved)
                BestXTable = this.XTrace(Iteration,:);
            else
                BestXTable = [];
                MinObserved = NaN;
                Iteration = NaN;
            end
        end
        
        function [BestXTable, MinMean] = minMeanPoint(this)
            % Attempts to find the feasible point BestXTable that mimimizes
            % the mean of the ObjectiveFcn model.
            if isempty(this.ObjectiveFcnGP) || ...
                    (this.PrivOptions.NumCoupledConstraints > 0 && any(cellfun(@isempty, this.ConstraintGPs)))
                BestXTable = [];
                MinMean = NaN;
            else
                % Set tighter optimization stopping criteria
                BOOptions = this.PrivOptions;
                BOOptions.MaxIterPerRestart = 40;
                BOOptions.RelTol = 1e-7;
                [MinMean, BestX] = findIncumbent(this, BOOptions);
                BestXTable = untransformPoints(this, BestX, true);
            end
            % If none found, revert to the minVisitedMeanPoint
            if isempty(BestXTable)
                [BestXTable, MinMean] = minVisitedMeanPoint(this);
            end
        end
        
        function [BestXTable, MinUCI] = minUCIPoint(this, Alpha)
            % Attempts to find the feasible point BestXTable that mimimizes
            % the 100(1-Alpha)% upper confidence interval (UCI) of the
            % ObjectiveFcn model.
            if isempty(this.ObjectiveFcnGP) || ...
                    (this.PrivOptions.NumCoupledConstraints > 0 && any(cellfun(@isempty, this.ConstraintGPs)))
                MinUCI = NaN;
                BestXTable = [];
            else
                BOOptions = this.PrivOptions;
                BOOptions.MaxIterPerRestart = 40;
                BOOptions.RelTol = 1e-7;
                VarSpec = BOOptions.VarSpec;
                BestX = iFminbndGlobal(@(X)GPF_UCI_OnPoints(this, X, Alpha), VarSpec.LBTrans, VarSpec.UBTrans, ...
                    BOOptions.NumRestartCandidates, BOOptions.NumRestarts, BOOptions.VerboseRestarts, ...
                    BOOptions.MaxIterPerRestart, BOOptions.RelTol);
                % XAtMin is not discretized (although it was discretized during optimization)
                BestX = legalizePoints(this, BestX);
                MinUCI = GPF_UCI_OnPoints(this, BestX, Alpha);
                BestXTable = untransformPoints(this, BestX, true);
            end
            % If none found, revert to the minVisitedUCIPoint
            if isempty(BestXTable)
                [BestXTable, MinUCI] = minVisitedUCIPoint(this, Alpha);
            end
        end
        
        function [BestXTable, MinMean, Iteration] = minVisitedMeanPoint(this)
            % Attempts to find the feasible point BestXTable that mimimizes
            % the mean of the ObjectiveFcn model, from among those points
            % visited (evaluated).
            if isempty(this.ObjectiveFcnGP) || ...
                    (this.PrivOptions.NumCoupledConstraints > 0 && any(cellfun(@isempty, this.ConstraintGPs)))
                BestXTable = [];
                MinMean = NaN;
                Iteration = NaN;
            else
                ObjMeans = GPFMeanOnPoints(this, this.XTrain);
                [MinMean, Iteration] = nanmin(ObjMeans);
                BestXTable = this.XTrace(Iteration,:);
            end
        end
        
        function [BestXTable, MinUCI, Iteration] = minVisitedUCIPoint(this, Alpha)
            % Attempts to find the feasible point BestXTable that mimimizes
            % the 100(1-Alpha)% upper confidence interval (UCI) of the
            % ObjectiveFcn model, from among those points visited
            % (evaluated).
            UCI = GPF_UCI_OnPoints(this, this.XTrain, Alpha);
            [MinUCI, Iteration] = nanmin(UCI);
            BestXTable = this.XTrace(Iteration,:);
        end
        
        function UCI = GPF_UCI_OnPoints(this, X, Alpha)
            % Return the Alpha% UCI of the ObjectiveFcnGP value for all
            % feasible rows of X. Return Inf for infeasible rows.
            UCI = inf(size(X,1),1);
            if isempty(this.ObjectiveFcnGP) || ...
                    (this.PrivOptions.NumCoupledConstraints > 0 && any(cellfun(@isempty, this.ConstraintGPs)))
                return;
            else
                [Xcanon, InputFeasible] = legalizePoints(this, X);
                allSat = allConstraintsSatisfied(this, Xcanon, this.ConstraintGPs, this.ErrorGP);
                feasible = InputFeasible & allSat;
                [FMean, YSD] = predict(this.ObjectiveFcnGP, Xcanon(feasible,:));
                FSD = funStd(YSD, this.ObjectiveFcnGP.Sigma);
                UCI(feasible) = norminv(1-Alpha, FMean, FSD);
            end
        end
        
        %% Misc
        function [X, this] = gridXFeasiblePoint(this)
            % A random grid point that is XFeasible, or [] if none exist.
            [X, this] = sampleGridWithoutReplacement(this);
            while ~isempty(X) && ~isInputFeasible(this, X)
                [X, this] = sampleGridWithoutReplacement(this);
            end
        end
        
        function X = randomXFeasiblePoints(this, N)
            % Generate N random initial points that satisfy the xConstraint. Return
            % then in transformed space.
            for i = 1:N
                X(i,:) = randomXFeasiblePoint(this);
            end
        end
        
        function XBest = randomXFeasiblePoint(this)
            X = iUniformPointsWithinBounds(this.PrivOptions.NumRestartCandidates, ...
                this.PrivOptions.VarSpec.LBTrans, this.PrivOptions.VarSpec.UBTrans);
            idx = find(isInputFeasible(this, X), 1);
            XBest = X(idx,:);
        end
        
        function XTable = conditionalizeX(this, X)
            XTable = untransformPoints(this, X, true);
            XTable = applyConditionalVariableFcn(this, XTable);
        end
        
        function NewXTable = applyConditionalVariableFcn(this, XTable)
            NewXTable = XTable;
            if ~isempty(this.PrivOptions.ConditionalVariableFcn)
                try
                    NewXTable = this.PrivOptions.ConditionalVariableFcn(XTable);
                catch me
                    bayesoptim.warn('ConditionalVariableFcnError');
                    rethrow(me);
                end
                % Check output
                if ~(istable(NewXTable) && ...
                     isequal(size(NewXTable), size(XTable)) && ...
                     isequal(varfun(@class, NewXTable, 'OutputFormat','cell'), varfun(@class, XTable, 'OutputFormat','cell')) && ...
                     isequal(NewXTable.Properties.VariableNames, XTable.Properties.VariableNames))
                    bayesoptim.err('ConditionalVariableFcnOutput');
                end
            end
        end
        
        function XTable = canonicalizePoints(this, XTable)
            % Map missing values to the minimum of their range or their first category.
            VarSpec = this.PrivOptions.VarSpec;
            missing = ismissing(XTable);
            for v = 1:width(XTable)
                if any(missing(:,v))
                    if VarSpec.isCat(v)
                        cats = VarSpec.Categories{v};
                        val = cats(1);
                    else
                        val = VarSpec.LBs(v);
                    end
                    XTable{missing(:,v), v} = val;
                end
            end
        end
        
        function X = transformPoints(this, XTable)
            % Convert a table in native space to a matrix of transformed points.
            VarSpec = this.PrivOptions.VarSpec;
            IsNumericLinear = ~VarSpec.isCat & ~VarSpec.isLog;
            IsNumericLog = ~VarSpec.isCat & VarSpec.isLog;
            if any(IsNumericLinear)
                X(:,IsNumericLinear) = XTable{:,IsNumericLinear};
            end
            if any(IsNumericLog)
                X(:,IsNumericLog) = log(XTable{:,IsNumericLog});
            end
            if any(VarSpec.isCat)
                X(:,VarSpec.isCat) = table2array(varfun(@intCodeCategorical, XTable(:,VarSpec.isCat)));
            end
        end
        
        function T = checkAndPrepareTableForPrediction(this, XTable, caller)
            % Check input, then find the required variables by Name and put
            % them in the order matching VarSpec.
            if isempty(XTable)
                T = table;
                return;
            elseif ~istable(XTable)
                bayesoptim.err('PredictArgNotTable', caller);
            else
                T = table;
                VarNames = this.PrivOptions.VarSpec.Names;
                for i = 1:numel(VarNames)
                    if ~ismember(VarNames{i}, XTable.Properties.VariableNames)
                        bayesoptim.err('PredictVarMissing', VarNames{i});
                    end
                    T.(VarNames{i}) = XTable.(VarNames{i});
                end
            end
        end
        
        function [KernelParams0, ConstantKernelParameters] = kernelParamsForFitrgp(this, TrainingData)
            VarSpec = this.PrivOptions.VarSpec;
            sigmaF = nanstd(TrainingData)/sqrt(2);
            if isnan(sigmaF) || sigmaF == 0
                sigmaF = 1;
            end
            KernelParams0 = [(VarSpec.UBTrans(:)-VarSpec.LBTrans(:))/2; sigmaF];    % Last one is the default kernel amplitude
            KernelParams0(VarSpec.isCat) = this.PrivOptions.CatKernelScale;
            ConstantKernelParameters = [VarSpec.isCat(:); false];
        end
    end
end

%% Local functions
function [VariableDescriptions, RemainingArgs] = iParseVariableDescriptionsArg(Args)
[VariableDescriptions,~,RemainingArgs] = internal.stats.parseArgs(...
    {'VariableDescriptions'}, {[]}, Args{:});
end

function val = iNanIfBad(val)
% val must be a finite real number. If not, set it to NaN.
if ~isscalar(val) || ~isnumeric(val) || ~isreal(val) || ~isfinite(val)
    val = NaN;
end
end

function GP = iFitrgpRobust(customOpts, X, Y, SigmaLowerBound, varargin)
% Fit a GP. If it fails due to singularity, iteratively double
% SigmaLowerBound until it succeeds (up to 10 times). If it never succeeds,
% return [];
if all(isnan(Y))
    bayesoptim.err('AllModelTargetsNaN');
end
SigmaLowerBound = max(SigmaLowerBound, 1e-6);
GP = [];
success = false;
doublings = 0;
C = bayesoptim.suppressWarnings();

if ischar(customOpts.covarianceFcn)
    kernelFcn = customOpts.covarianceFcn;
else
    kernelFcn = @(Xm, Xn, theta) customOpts.covarianceFcn(Xm, Xn, theta, customOpts);
end

for i=1:size(varargin,2)
    if strcmp('KernelParameters', varargin{i})
        kParams = varargin{i+1}(end-1:end);
    end
end

varargin{end+1} = 'ConstantKernelParameters';
varargin{end+1} = logical([0;0]);
varargin{end+1} = 'KernelParameters';
varargin{end+1} = kParams;
varargin{end+1} = 'KernelFunction';
varargin{end+1} = kernelFcn;
% if size(Y,2) > 3
%     varargin{end+1} = 'OptimizeHyperparameters';
%     varargin{end+1} = {'KernelScale'};
%     varargin{end+1} = 'HyperparameterOptimizationOptions';
%     varargin{end+1} = struct('ShowPlots',1,'Verbose',1); %, 'UseParallel', 1);
% end
% varargin{end+1} = 'Standardize';
% varargin{end+1} = true;

while ~success && doublings <= 10
    try
        GP = compact(fitrgp(X, Y, varargin{:}, 'SigmaLowerBound', SigmaLowerBound));
        success = true;
    catch me
        SigmaLowerBound = 2*SigmaLowerBound;
        doublings = doublings + 1;
    end
end
if ~success
    xxx = 1;
end
end

function s = iGetSigmaF(GPR)
s = GPR.KernelInformation.KernelParameters(end);
end

function [x, fval] = iFminbndGlobal(fcn, LB, UB, NumCandidates, BestK, ...
    Verbose, MaxIter, TolFun)
% Attempts to find the global minimum of fcn() within bounds given by the
% vectors LB and UB. Generates NumCandidates random points within bounds,
% runs fcn() on all of them and selects the BestK best. Then runs a bounded
% fminsearch on each of the K points and returns the best overall. The
% flags Verbose, MaxIter, and TolFun are options passed directly to
% fminsearch.

% Choose the best K points from N random candidates
[x0, fvals] = iMinKPointsFromN(fcn, LB, UB, NumCandidates, BestK);
% Do a local optimization on each of the K points
% if sum(isfinite(fvals)) < BestK
%     bayesoptim.err('iMinKPointsFromNFailed');
% end
if Verbose
    VerboseOpts = {'Display', 'iter'};
else
    VerboseOpts = {'Display', 'off'};
end
for row = BestK:-1:1
    [xs(row,:), fvals(row)] = fminsearch(@boundedFcn, x0(row,:), ...
        optimset(...
        'MaxIter', MaxIter, ...
        'TolX', Inf, ...
        'TolFun', TolFun, ...
        VerboseOpts{:}));
end
[fval, loc] = min(fvals);
x = xs(loc,:);

    function y = boundedFcn(x)
        if any(x < LB | x > UB)
            y = Inf;
        else
            y = fcn(x);
        end
    end
end

function [X, fvals] = iMinKPointsFromN(fcn, LB, UB, N, K)
% Generate a random set of N points within bounds. Return the K points that
% best minimize fcn.
X = iUniformPointsWithinBounds(N, LB, UB);
fvals = fcn(X);
[~, rows] = sort(fvals, 'ascend');
X = X(rows(1:K), :);
fvals = fvals(rows(1:K));
end

function X = iUniformPointsWithinBounds(N, LB, UB)
Unifs = rand(N, size(LB,2));
X = Unifs .* repmat(UB - LB, N, 1);
X = X + repmat(LB, N, 1);
end

function tf = iAnyInitializationArgs(NVPs)
Names = {'InitialX', 'InitialObjective', 'InitialConstraintViolations', ...
    'InitialErrorValues', 'InitialObjectiveEvaluationTimes', 'InitialIterationTimes', ...
    'InitialUserData', 'NumSeedPoints'};
Defaults = cell(1, numel(Names));
[~,~,~,~,~,~,~,~,setflag,~] = internal.stats.parseArgs(Names, Defaults, NVPs{:});
tf = bayesoptim.anyFlagsSet(setflag);
end

function Indices = iFindUnusedGridIndices(NumGridDivisions, SampledGridIndices)
if isempty(SampledGridIndices)
    for p=numel(NumGridDivisions):-1:1
        Indices(p) = randi(NumGridDivisions(p));
    end
else
    Indices = SampledGridIndices(1,:);
    while ismember(Indices, SampledGridIndices, 'rows')
        for p=1:numel(NumGridDivisions)
            Indices(p) = randi(NumGridDivisions(p));
        end
    end
end
end

function checkVariableDescriptions(NewVariableDescrips, OldVariableDescrips)
% The variable names to be optimized must be the same. Numeric variables
% may have their Range, Type, and Transform changed, but must remain
% numeric. Categoricals cannot be changed.
Old = sortVDByName(OldVariableDescrips([OldVariableDescrips.Optimize]));
New = sortVDByName(NewVariableDescrips([NewVariableDescrips.Optimize]));
if ~isequal({Old.Name}, {New.Name})
    bayesoptim.err('ResumeVarsChanged');
end
for i = 1:numel(Old)
    switch Old(i).Type
        case 'categorical'
            if ~isequal(New(i).Type, 'categorical') || ~isequal(Old(i).Range, New(i).Range)
                bayesoptim.err('ResumeCatsChanged');
            end
        case {'integer','real'}
            if ~ismember(New(i).Type, {'integer','real'})
                bayesoptim.err('ResumeNumericChanged');
            end
    end
end
end

function VD = sortVDByName(VD)
[~,I] = sort({VD.Name});
VD = VD(I);
end

function doubleVec = intCodeCategorical(catVec)
doubleVec = double(catVec);
end

function ObjectiveNargout = nargoutFromNumConstraints(NumCoupledConstraints)
% We know nargout < 3
if NumCoupledConstraints > 0
    ObjectiveNargout = 2;
else
    ObjectiveNargout = 1;
end
end

function S = nancumsum(X)
X(isnan(X)) = 0;
S = cumsum(X);
end

function FSD = funStd(YSD, Sigma)
% Return the standard deviation of the posterior over function values (F),
% given the std of the posterior over observed values (Y) and the noise
% standard deviation (Sigma). FSD, YSD, and Sigma are arrays.
FSD = sqrt(max(0, YSD.^2 - Sigma.^2)); % Use max to avoid complex numbers.
end

function Cleanup = holdOn(Axes)
% Turn hold on and return an oncleanup object to turn hold off.
hold(Axes, 'on');
Cleanup = onCleanup(@()holdOff(Axes));
    function holdOff(A)
        if isvalid(A)
            hold(A, 'off');
        end
    end
end

