function [A,b,vert] = voronoiPolyhedrons( seeds,varargin )
%Compute inequalities and, if desired, the vertices of Euclidean voronoi
%polyhedrons.
%
%    [A,b,vert] = voronoiPolyhedrons( seeds,lb,ub )
%          
%IN:        
%           
%    seeds: an NxM matrix of N-dimensional seed points
%       lb: N-vector of lower bounds
%       ub: N-vector of upper bounds
%           
%OUT:       
%           
%       A,b: length M cell arrays such that A{i}*x<=b{i} are the inequalities  
%            of the i-th Voronoi polyhedron (i.e., of seed(:,i)), but
%            restricted to bounds lb,ub
%
%    vert: length M cell array such that vert{i} are the vertices of the
%          i-th polyhedron. For this output to be requested, finite box bounds, lb,
%          ub must be given.
%
% source 
% https://www.mathworks.com/matlabcentral/answers/258667-vertices-of-polygons-in-bounded-voronoi-diagram
%


  seeds=seeds.';
  [N,d]=size(seeds);

  psquares=sum(seeds.^2,2);
  bgraph=bsxfun(@minus, psquares(:), psquares(:).')/2;
  
  Agraph=bsxfun(@minus,seeds,reshape(seeds.',1,d,[]));
  
  
  A=num2cell(Agraph,[1,2]); A=A(:).';
  b=num2cell(bgraph,1); b=b(:).';
  
  
  if nargin>1 %Include bounds
      
    [Aq,bq]=addBounds([],[],[],[],varargin{:});

  else
      
      Aq=[]; bq=[];
    
  end 
  
  
  
      for i=1:N

          A{i}(i,:)=[];
          b{i}(i)=[];
          
          A{i}=[A{i};Aq];
          b{i}=[b{i};bq];

      end



  
  if nargout>2 %Compute vertices
      
      vert=cell(1,N);
     
      if nargin==1, error 'Bounds must be supplied for vertex computation'; end
      
       for i=1:N

          vert{i}=qlcon2vert(seeds(i,:),A{i},b{i}).';

       end   
      
  end
