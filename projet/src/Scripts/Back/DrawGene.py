import matplotlib.pyplot as plt

"""
MIT License
Copyright (c) [2016] [Parashar Dhapola]
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

__author__ = "Parashar Dhapola"
__email__ = "parashar.dhapola@gmail.com"


class GeneImage(object):
    def __init__(self, exon_intervals, marker_pos=[], marker_heights=[], marker_colors=['red', 'blue'],
                 marker_size=100, marker_weight=1.5, exon_color="black", intron_color="grey",
                 intron_weight=2, intron_style='-', bar_color='cornflowerblue', bar_xmin : int = None, bar_xmax : int = None,
                 bg_color="white", show_labels=True):
        self.canvas : plt.Axes = None
        self.exonIntervals = exon_intervals
        self.markerPositions = marker_pos
        self.markerHeights = marker_heights
        self.markerColors = marker_colors
        self.markerSize = marker_size
        self.MarkerWeight = marker_weight
        self.exonColor = exon_color
        self.intronColor = intron_color
        self.intronWeight = intron_weight
        self.intronStyle = intron_style
        self.barColor= bar_color
        self.barColorXmin = self.exonIntervals[0][0] if bar_xmin is None else bar_xmin
        self.barColorXmax = self.exonIntervals[-1][1] if bar_xmax is None else bar_xmax
        self._check_input()
        self.bgColor = bg_color
        self.markerDefaultColor = 'grey'
        self.numExons = len(self.exonIntervals)
        self.totalSpan = self.exonIntervals[-1][1] - self.exonIntervals[0][0]
        self.minExonLen = self.totalSpan*0.005
        self.ylims = {'exon_max': 2, 'exon_min':1}
        self.figure, self.canvas = plt.subplots(figsize=(15,1.5))
        self.show_labels = show_labels
        self._draw()

    def _check_input(self):
        if not isinstance(self.barColorXmin, int) or not isinstance(self.barColorXmax, int):
            raise ValueError("Border values must be integers")
        if not isinstance(self.markerPositions, list) or not isinstance(self.markerColors, list):
            raise ValueError("Markers positions and colors must be lists")
        if not isinstance(self.exonIntervals, list):
            raise ValueError("Exon intervals must be a list")
        for intervals in self.exonIntervals:
            if not isinstance(intervals, tuple) or len(intervals) != 2:
                raise ValueError("Exon intervals must be a list of lists with 2 values")
            if not isinstance(intervals[0], int) or not isinstance(intervals[1], int):
                raise ValueError("Exon intervals values must be integers")

    def _set_limits(self):
        self.ylims['intron_max'] = self.ylims['exon_max']*0.9
        self.ylims['intron_min'] = (self.ylims['exon_max'] + self.ylims['exon_min'])/2.0
        self.ylims['bar_min'] = self.ylims['exon_max']+0.2
        self.ylims['bar_max'] = self.ylims['bar_min']+(self.ylims['exon_max']-self.ylims['exon_min'])/5.0
        
    
    def _transform_spans(self):
        span_lens = [x[1]-x[0] for x in self.exonIntervals]
        max_len = float(max(span_lens))
        transformed_intervals = []
        if max_len < self.minExonLen:
            span_ratios = [x/max_len for x in span_lens]
            expansion_factor = self.totalSpan*1e-11
            for i in range(1,10):
                ef = (2**i)*expansion_factor
                if max_len+ef > self.minExonLen:
                    expansion_factor = ef
                    break
            for i,j in zip(self.exonIntervals, span_ratios):
                mid = (i[0] + i[1])/2
                f = (expansion_factor*j)/2
                if mid+f - mid-f > self.minExonLen:
                    transformed_intervals.append([mid-f, mid+f])
                else:
                    transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
        else:
            for i in range(self.numExons):
                if span_lens[i] < self.minExonLen:
                    mid = (self.exonIntervals[i][0] + self.exonIntervals[i][0])/2 
                    transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
                else:
                    transformed_intervals.append(self.exonIntervals[i])
        self.exonIntervals = transformed_intervals[:]
        
    def _draw_exon(self, span):
        self.canvas.fill_between(span, self.ylims['exon_min'], self.ylims['exon_max'],
                                 edgecolor=self.bgColor, facecolor=self.exonColor)
        return True
        
    def _draw_intron(self, span):
        mid = (span[0]+span[1])/2.0
        self.canvas.plot([span[0], mid], [self.ylims['intron_min'], self.ylims['intron_max']],
                         c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
        self.canvas.plot([mid, span[1]], [self.ylims['intron_max'], self.ylims['intron_min']],
                         c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
        return True
    
    def _draw_markers(self):
        if self.markerHeights == []:
            self.markerHeights = [self.ylims['exon_max']-self.ylims['exon_min'] for x in self.markerPositions]
        if self.markerColors == []:
            self.markerColors = [self.markerDefaultColor for x in self.markerPositions]           
        for p,h,c in zip(self.markerPositions, self.markerHeights, self.markerColors):
            self.canvas.plot((p, p), (self.ylims['bar_max'], self.ylims['bar_max']+h),
                             linestyle='-', color='black', linewidth=self.MarkerWeight, alpha=0.7)
            self.canvas.scatter(p, self.ylims['bar_max']+h+0.25, s=self.markerSize, marker='o', c=c,
                                edgecolor=c, alpha=1)
        
    
    def _clean_axes(self):
        #self.canvas.set_ylim((self.ylims['exon_min'], self.ylims['bar_max']))
        self.canvas.set_yticks([], [])
        self.canvas.get_xaxis().tick_top()
        self.canvas.tick_params(axis='x', direction='out')
        self.canvas.set_xticks([])
        for o in ["top", "bottom", "left", "right"]:
            self.canvas.spines[o].set_visible(False)
        min_pos = int(self.exonIntervals[0][0] - self.totalSpan * 0.1)
        if min_pos < 0:
            min_pos = 0
        max_pos = int(self.exonIntervals[-1][1] + self.totalSpan * 0.1)
        minortick_pos = [x for x in range(min_pos, max_pos, (max_pos-min_pos)//20)][1:]
        for i in minortick_pos:
            self.canvas.axvline(i, alpha=0.2, c='black', ls='--')
        self.canvas.text(minortick_pos[0], self.ylims['exon_min']-0.5,
                         minortick_pos[0], fontsize=8, ha='center')
        self.canvas.text(minortick_pos[-1], self.ylims['exon_min']-0.5,
                         minortick_pos[-1], fontsize=8, ha='center')
        self.canvas.set_xlim(minortick_pos[0]-(minortick_pos[1]-minortick_pos[0]),
                             minortick_pos[-1]+(minortick_pos[-1]-minortick_pos[-2]))
        
    def _add_labels(self):
        """
        Method that add the distances below the blue square
        """
        dist = self.barColorXmax - self.barColorXmin
        self.canvas.text(x = self.barColorXmin + dist//2, y = self.ylims['bar_min']+0.5, s = str(dist), fontsize=8, ha='center')
        
    def _draw(self):
        self._set_limits()
        self._transform_spans()
        for i in range(self.numExons):
            if i > 0:
                self._draw_intron([self.exonIntervals[i-1][1], self.exonIntervals[i][0]])
            self._draw_exon(self.exonIntervals[i])
        self.canvas.fill_between([self.barColorXmin, self.barColorXmax],
                                  self.ylims['bar_min'], self.ylims['bar_max'],
                                  edgecolor=self.bgColor, facecolor=self.barColor)
        self._draw_markers()
        self._clean_axes()
        if self.show_labels:
            self._add_labels()
    
    def show(self):
        plt.show()
        # plt.close()
    
    def save(self, path):
        plt.savefig(path, bbox_inches='tight', pad_inches=0.1)
        plt.close()
        
if  __name__ == "__main__":
    #exons positions
    exon_pos = [
        [97543299, 97544702], [97547885, 97548026], [97564044, 97564188],
        [97658624, 97658804], [97700407, 97700550], [97770814, 97770934]
    ]
    #marker positions
    marker_pos = [97647885, 98247485]
    gene = GeneImage(exon_pos, marker_pos, marker_colors=['red', 'blue'], bar_xmin=97647885, bar_xmax=98247485 )
    gene.show()