classdef ISINGrid
    % ISIN grid
    
    
    % Class properties
    properties
        rowCount
        lats
        deltaLat
        rowLength
        deltaLons
        binOffsets
        totalBinCount
        mapX
        mapY
    end  %properties
    
    % Class methods
    methods
        function obj = ISINGrid(rc)
            if nargin > 0
                obj.rowCount = rc;
            else
                obj.rowCount = 4320;  % 4 km grid
            end
            
            obj.deltaLat = 180.0 / obj.rowCount;
            
            %         /* effective size (km and deg) of a cell in latitudinal direction */
            %         /* computes the number of cells in longitudinal direction for each rowIndex */
            %         /* and the bin effective longitudinal size */
            %         lats = new double [rowCount];
            %         deltaLons = new double [rowCount];
            %         rowLength = new int [rowCount];
            %         binOffsets = new int [rowCount];
            binCount = 0;
            for i = 1:obj.rowCount
                j = i - 1;
                obj.lats(i) = -90.0 + (i - 0.5) * obj.deltaLat;
                obj.rowLength(i) =  fix(0.5 + 2.0 * obj.rowCount * cos(obj.lats(i)*pi/180.));
                obj.deltaLons(i) = 360.0 / obj.rowLength(i);
                if (i == 1)
                    obj.binOffsets(i) = 0;
                else
                    obj.binOffsets(i) = obj.binOffsets(i - 1) + obj.rowLength(i - 1);
                end
                binCount = binCount + obj.rowLength(i);
            end
            obj.totalBinCount = binCount;
        end  % constructor
        
        function binOffset =  getBinOffset(obj, rowIndex)
            binOffset = obj.binOffsets(rowIndex);
        end
        
        function rowIndex = getRowIndex(obj, lat)
            rowIndex = fix((lat + 90.0)*(obj.rowCount)/180)+1;
            rowIndex = min(rowIndex, obj.rowCount);
        end
        
        function rowIndex = getRowIndex2(obj, binIndex)
            g = find(binIndex >= obj.binOffsets);
            rowIndex = g(end);
        end
        
        function rowLength = getRowLength(obj, rowIndex)
            rowLength = obj.rowLength(rowIndex);
        end
        
        %     /**
        %      * Gets the zero-based column index in the global ISIN grid for the given zero-based row index and longitude.
        %      *
        %      * @param rowIndex the zero-based row index in the range 0...{@link #getRowCount()}-1
        %      * @param lon      the longitude in the range 0...360 degree or -180...+180
        %      *
        %      * @return the zero-based column index only if the longitude is in the range 0...360 degree, otherwise undefined
        %      *
        %      */
        function colIndex = getColIndex(obj, rowIndex, lon)
            if (lon >= 180)
                lon = lon - 360;
            end
            colIndex =  fix((lon + 180.0) * obj.getRowLength(rowIndex)/360.0)+1;
        end
        
        function binIndex = getBinIndexLon(obj, rowIndex, lon)
            colIndex = obj.getColIndex(rowIndex, lon);
            binIndex = obj.getBinOffset(rowIndex) + colIndex;
        end
        
        function binIndex = getBinIndex(obj, lat, lon)
            rowIndex = obj.getRowIndex(lat);
            binIndex = obj.getBinIndexLon(rowIndex, lon);
        end
        
        
        function lat = getLat(obj, rowIndex)
            lat = obj.lats(rowIndex);
            return;
        end

        function p = getGridPoint(obj, binIndex)
            rowIndex = obj.getRowIndex2(binIndex);
            if rowIndex == -1
                colIndex = -1;
            else
                colIndex = binIndex - obj.binOffsets(rowIndex);
            end
            p = zeros(2,1);
            p(1) = obj.getLat(rowIndex);
            p(2) = ((360.0*(colIndex + 0.5)/obj.rowLength(rowIndex))-180.0);
            return;
        end
        
        function h = map(obj, lats, lons, chl)
            figure;
            m_proj('sinusoidal');
            %             m_gshhs_('save','world');
            if (isempty(obj.mapX))
                [obj.mapX, obj.mapY] = m_ll2xy(lons, lats);
            end
            m_proj('sinusoidal', 'lon', [-180 180], 'lat', [-90 90], 'clongitude', 0);
            h = color_mark(obj.mapX, obj.mapY, log10(chl));
            set(h, 'markersize', 3);
            m_grid;
            m_coast('line', 'color','k');
            %             m_grid('xtick',[-90 -60 -30],'tickdir','out','ytick',[30 60],'linest','-','color', linecol, 'backcolor', backcol);
            %             m_usercoast('azmpl','patch', fillcol, 'edgecolor', 'none');
            
        end
        
    end  %methods
end  % classdef
