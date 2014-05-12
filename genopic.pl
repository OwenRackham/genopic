#!/usr/bin/env perl

use warnings;
use strict;
use POSIX;
use Math::Complex;
=head1 NAME

B<genopic.cgi> - Display the architecture for a specified SUPERFAMILY protein.

=head1 DESCRIPTION

Outputs an SVG rendering of the given proteins structual and disordered architecture. Weaker hits are included with their e-values specified as 'hanging' blocks.

An example use of this script is as follows:

To emulate SUPERFAMILY genome page style figures as closely as possible include something similar to the following in the page:

<div width="100%" style="overflow:scroll;">
	<object width="100%" height="100%" data="/cgi-bin/disorder.cgi?proteins=3385949&genome=at&supfam=1&ruler=0" type="image/svg+xml"></object>
</div>

To have super duper Matt style figures do something like:

<div width="100%" style="overflow:scroll;">
	<object width="100%" height="100%" data="/cgi-bin/disorder.cgi?proteins=3385949,26711867&callouts=1&ruler=1&disorder=1" type="image/svg+xml"></object>
</div>


=head1 TODO

B<HANDLE PARTIAL HITS!>

I<SANITIZE INPUT MORE!>

	* Specify lists of proteins, along with other search terms like comb string, required by SUPERFAMILY.

=head1 AUTHOR

B<Owen Rackham> - I<owen.rackham@bristol.ac.uk>

=head1 NOTICE

B<Owen Rackham> (2014) First features added.

=head1 LICENSE AND COPYRIGHT

B<Copyright 2014 Owen Rackham>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 FUNCTIONS

=over 4

=cut

use POSIX qw/ceil floor/;
use Data::Dumper;
use JSON::XS;

#Global variables
my $filename = $ARGV[0];
my $output = '';
my $scale = 1;
my $forprint = 0;
my $force_png = 0;
my $xpad = 10;
my $ypad = 10;

my ($genos,$count) = read_23andMe($filename);
plot_genotypes($genos,$count);
=item B<read_23andMe>
    read the 23andMe and file and store genotypes in an array
=cut
sub read_23andMe {
    my $file = shift;
    my @genotypes;
    my $count = 0;
    open FILE,"$file";
    while(<FILE>){
        my $line = $_;
        unless($line =~ '#'){
            my @line = split("\t",$line);
            push(@genotypes,$line[3]);
            $count++;
        }
    }
    return(\@genotypes,$count);
}

=item B<hsv2rgb($Hue, $Saturation, $Value)>

	Function to convert HSV colour space values to RGB colour space.
	Returns RGB value as [R,G,B]
=cut
sub hsv2rgb {
	my ($Hue,$Saturation,$Value) = @_;
	my ($Red,$Green,$Blue) = (0,0,0);
	
	#Check the input and warn if it's a bit wrong
	warn "Invalid Hue component of HSV colour passed, with value: $Hue." unless ($Hue >= 0.0 and $Hue <= 360.0);
	warn "Invalid Saturation component of HSV colour passed, with value: $Saturation." unless($Saturation >= 0.0 and $Saturation <= 1.0);
	warn "Invalid Value component of HSV colour passed, with value: $Value." unless ($Value >= 0.0 and $Value <= 1.0);
	
	#If colour has no saturation return greyscale RGB
	if ($Saturation == 0) {
		$Red = $Green = $Blue = $Value;
		return [$Red, $Green, $Blue];
	}
	
	#Partition the Hue into the 5 different colour chroma and then map each of these to RGB based on the colour theory
	$Hue /= 60.0;
	my $Chroma = floor($Hue) % 6; 
	my $H_d = $Hue - $Chroma; 
	
	#RGB cube components
	my ($I,$J,$K) = ( $Value * ( 1 - $Saturation ),
	                            $Value * ( 1 - $Saturation * $H_d ),
				    $Value * ( 1 - $Saturation * ( 1 - $H_d ) )
                                    );
	
	#Map components to RGB values per chroma
	if ($Chroma == 0) { ($Red,$Green,$Blue) = ($Value,$K,$I); }
	elsif ($Chroma == 1) { ($Red,$Green,$Blue) = ($J,$Value,$I); }
	elsif ($Chroma == 2) { ($Red,$Green,$Blue) = ($I,$Value,$K); }
	elsif ($Chroma == 3) { ($Red,$Green,$Blue) = ($I,$J,$Value); }
	elsif ($Chroma == 4) { ($Red,$Green,$Blue) = ($K,$I,$Value); }
	else{ ($Red,$Green,$Blue) = ($Value,$I,$J); }
	
	#Return the RGB value in the integer range [0,255] rather than real [0,1]
	return [floor($Red * 255),floor($Green * 255),floor($Blue * 255)];
}

#Internal sub to randomly shuffle an array
sub fisher_yates_shuffle {
    my ($array) = @_;
    my $current;
    for ($current = @$array; --$current; ) {
        my $selected = int rand ($current+1);
        next if $current == $selected;
        #Reverse the slice between current position and the randomly selected
        @$array[$current,$selected] = @$array[$selected,$current];
    }
    return $array;
}

=item B<colourset($num_colours, $method)>

	Function to grab a set of well spaced colours from an HSV colour wheel.
		$num_colours - The number of colour values to produce, must be greater than 0 but no bigger than 360
		$method - The method for selecting colours over HSV colour space, either 'equal_spacing' or for around 10 colours 'chroma_bisection' is better.
	Returns an array of RGB values of the form ([R,G,B]) and undef on $num_colours out of bounds
=cut
sub colourset {
	my ($num_colours,$method,$value) = @_;
	if ($num_colours <= 0 or $num_colours > 360) {
		warn "Number of colours requested out of bounds.";
		return undef;
	}
	$method = 'chroma_bisection' unless $method;
	
	
	#Colours to return
	my %colours;
	
	#Default Hue Saturation and Value, saturation of 0.65 gives a more pastel feel!
	my ($Hue, $Saturation, $Value) = (0.0,0.65,0.95);
	$Value = $value if defined $value;

	#The interval to space colours around the wheel if equal
	my $hsv_interval = 360 / $num_colours;
	
	#Array of degrees for reuse to create ranged arrays with a given interval
	my @degrees = 1..360;
	
	#Iteratively bisect each chroma segment so that the first 6 colours are well spaced perceptually.
	#However after 12 colours we will have increasing pairs that are more confused as 
	#they are increasingly close to each other compared to the rest of the colours!
	#To get around this problem of selecting closely around a single bisection, we jump around the 
	#chroma randomly sampling.
	if ($method eq 'chroma_bisection') {
		#The current cycle of chroma bisection
		my $hsv_cycle = 1;
		#Number of colours selected by bisecting chroma so far
		my $colours_selected = 0;
		
		#While we still have colours to select
		while ($colours_selected != $num_colours) {
			#Work out the size of interval to use this cycle around the wheel
			$hsv_interval = 60 / $hsv_cycle;
			#Get all the Hues for this cycle that haven't already been examined and are on the line of bisection
			my @Hues = grep { not exists $colours{$_%360} } range({ from => 1, to => 360, by => $hsv_interval});
			#Shuffle so that we don't take from around the same chroma all the time, only perceptually worthwhile after 12th colour
			fisher_yates_shuffle(\@Hues) if $hsv_cycle > 2;
			
			#While we still have hues to select from in this cycle
			while (@Hues) {
				#Finish if we have enough colours
				last if $colours_selected == $num_colours;
				#Consume a Hue from this cycle
				$Hue = shift @Hues;
				#360 should be 0 for red
				$Hue %= 360;
				#Store this Hue and mark selection
				$colours{$Hue} = hsv2rgb($Hue,$Saturation,$Value) ;
				$colours_selected++;
			}
			$hsv_cycle++;
		}
	}
	
	#Just space colours even distances apart over the HSV colour wheel.
	#You have slightly odd/garish colours coming out, but you dont get uneven perceptual distance
	#between pairs of colours. This scales far better despite the horrible colours.
	elsif ($method eq 'equal_spacing') {	
		foreach $Hue (1..$num_colours) {
			$Hue = ($Hue * $hsv_interval) % 360;
			$colours{$Hue} = hsv2rgb($Hue,$Saturation,$Value) ;
		}
	}
	
	#Otherwise return nothing and warn the programmer
	else {
		warn "Colourset method not known, use either 'equal_spacing' or for fewer colours 'chroma_bisection'";
		return undef;
	}
	
	#Shuffle final colours so that even if we do use chroma_bisection closer colours will hopefully not be sequential
	@_ = values %colours;
	#fisher_yates_shuffle(\@_);
	return @_;
}


=item B<header>

	Output the SVG header with some JavaScript for performing mouse popup information
=cut
sub header {
	my ($width, $height) = @_;
    $width *= $scale;
    $height *= $scale;
	my $header = <<EOF;
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
	 "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     width="$width" height="$height"
     viewBox="0 0 $width $height"
     preserveAspectRatio="xMidYMid meet"
     onload="init(evt)">
EOF

    unless ($forprint or $force_png) {
    $header .= <<EOF;
     <script type="text/ecmascript"><![CDATA[
      var svg_document, svg_root;
      var tooltop, tip_box, tip_title, tip_text;
      var hit_range, hit_quality;

      function init(evt) {
	 svg_document = evt.target.ownerDocument;
	 svg_root = document.documentElement;

	 tooltip = svg_document.getElementById('tooltip');
	 tip_box = svg_document.getElementById('tip_box');
	 tip_text = svg_document.getElementById('tip_text');
	 tip_title = svg_document.getElementById('tip_title');
	 tip_desc = svg_document.getElementById('tip_desc');
	 hit_range = svg_document.getElementById('hit_range');
	 hit_quality = svg_document.getElementById('hit_quality');

      };

      function show_tip(evt) {
	  var x = evt.clientX + window.pageXOffset;
	  var y = evt.clientY + window.pageYOffset;
	  var hit = evt.target;

	  var title_text = hit.getElementsByTagName('name').item(0);
	  if (! title_text ) {
	     return false;
	  }

	  title_text = title_text.firstChild.nodeValue;
	  tip_title.firstChild.nodeValue = title_text;
	  tip_title.setAttributeNS(null, 'display', 'inline' );

	  var description_text = hit.getElementsByTagName('desc').item(0);
	  if (description_text) {
	     description_text = description_text.firstChild.nodeValue;
	     tip_desc.firstChild.nodeValue = description_text;
	     tip_desc.setAttributeNS(null, 'display', 'inline' );
	  }
	  else {
	     tip_desc.setAttributeNS(null, 'display', 'none' );
	  }

	  var range_text = hit.getElementsByTagName('range').item(0);
	  if (range_text) {
	     range_text = range_text.firstChild.nodeValue;
	     hit_range.firstChild.nodeValue = range_text;
	     hit_range.setAttributeNS(null, 'display', 'inline' );
	  }
	  else {
	     hit_range.setAttributeNS(null, 'display', 'none' );
	  }

	  var quality_text = hit.getElementsByTagName('quality').item(0);
	  if (quality_text) {
	     quality_text = quality_text.firstChild.nodeValue;
	     hit_quality.firstChild.nodeValue = quality_text;
	     hit_quality.setAttributeNS(null, 'display', 'inline' );
	  }
	  else {
	     hit_quality.setAttributeNS(null, 'display', 'none' );
	  }

	  tip_title.firstChild.nodeValue = title_text;
	  tip_desc.firstChild.nodeValue = description_text;
	  hit_range.firstChild.nodeValue = range_text;
	  hit_quality.firstChild.nodeValue = quality_text;
	  
	  var box = tip_text.getBBox();
	  tip_box.setAttributeNS(null, 'width', Number(box.width) + 10);
	  tip_box.setAttributeNS(null, 'height', Number(box.height) + 10);
	  
	  tooltip.setAttributeNS(null, 'transform', 'translate(' + x + ',' + y + ')');
	  tooltip.setAttributeNS(null, 'visibility', 'visible');
      };


      function hide_tip(evt) {
	  tooltip.setAttributeNS(null, 'visibility', 'hidden');
      };

   ]]></script>
EOF
    }

    $header .= <<EOF;
  <defs>
    <linearGradient id="disorder"
		    x1="0%" y1="0%"
		    x2="0%" y2="100%"
		    spreadMethod="pad">
      <stop offset="0%"   stop-color="#cccccc" stop-opacity="0.6"/>
      <stop offset="100%" stop-color="#666666" stop-opacity="0.6"/>
    </linearGradient>
    <radialGradient id="radial-glow"
            fx="40%"
            fy="40%"
            r="55%"
            spreadMethod="pad">
        <stop offset="0%"   stop-color="#cccccc" stop-opacity="0.5" />
        <stop offset="100%" stop-color="#cccccc" stop-opacity="0.01" />
    </radialGradient>
    <filter id="emboss" >
        <feGaussianBlur in="SourceAlpha" stdDeviation="2" result="blur"/>
        <feSpecularLighting in="blur" surfaceScale="-3" style="lighting-color:white" specularConstant="1" specularExponent="16" result="spec" kernelUnitLength="1" >
            <feDistantLight azimuth="45" elevation="45" />
        </feSpecularLighting>
        <feComposite in="spec" in2="SourceGraphic" operator="in" result="specOut"/>
    </filter>
        <pattern id="binding" x="0" y="0" width="5" height="5"
        patternUnits="userSpaceOnUse">
        <path d="M 0 0 Q .25 5 2.5 2.5 T 5 5"
            style="stroke: black; fill: none;"/>
        </pattern>
  </defs>
  <g transform="translate($xpad,$ypad) scale($scale)">
EOF
}

=item B<footer>

	Print out the SVG footer containing the tooltip popup XML
=cut
sub footer {
	my $footer = <<EOF;
   </g>
EOF
	unless ($force_png or $forprint) {
	$footer .= <<EOF;
   <g id='tooltip' opacity='0.8' visibility='hidden' pointer-events='none'>
      <rect id='tip_box' x='0' y='5' width='88' height='20' rx='2' ry='2' fill='white' stroke='black'/>
      <text id='tip_text' x='5' y='20' font-family='Arial' font-size='10'>
        <tspan id='tip_title' x='5' font-weight='bold' text-decoration="underline"><![CDATA[]]></tspan>
	    <tspan id='hit_range' x='5' dy='15' fill='red'><![CDATA[]]></tspan>
        <tspan id='hit_quality' x='5' dy='10' fill='green'><![CDATA[]]></tspan>
	    <tspan id='tip_desc' x='5' dy='10' fill='blue'><![CDATA[]]></tspan>
      </text>
   </g>
EOF
    }
	$footer .= <<EOF;
</svg>
EOF
	return $footer;
}

sub plot_genotypes {
    my ($genotypes) = @_;
    my $genos = scalar(@{$genotypes});
    my $x = 0;
    my $y = 0;

    #Print out the SVG footer containing the tooltip popup markup/design

    #Print out SVG and ECMA script popup header
    my $width = 5940;
    my $height = 8410;
    my $squares = $width*$height;
    my $squares_per_gen = $squares/$genos;
    my $genosize = floor(sqrt($squares_per_gen));
    my $output;
    foreach my $gen (@{$genotypes}){
        $output .= "<rect x=\"$x\" y=\"$y\" width=\"$genosize\" height=\"$genosize\" style=\"fill:rgb(0,0,255);\"/>";
        if(($x+$genosize)<$width){
            $x = $x+$genosize;
        }else{
            $y = $y+$genosize;
            $x = 0;
        }
    }
    $output .= footer();    
    $output = header($width,$height) . $output;
    

    if ($force_png) {
	    use File::Temp qw/tempfile/;
	    my ($fh,$filename) = tempfile( 'architectureXXXXXX', UNLINK => 1, TMPDIR => 1, SUFFIX => '.svg' );
	    print $fh $output;
        #Using ImageMagick to convert to PNG
	    #open (PNG, '|-', "convert -background None -size ${width}x${height} svg:$filename --compress None --depth 32 png:fd:1") or error("Couldn't convert to PNG.");

        #Using RSVGlib directly, much nicer results!
        #If we want print quality up the DPI and pixel content
        if ($forprint) {
            $width *= $scale * 2;
            $height *= $scale * 2;
	        open (PNG, '|-', "rsvg -d 300 -p 300 --format=png -w ${width} -h ${height} $filename /dev/stdout") or error("Couldn't convert to PNG.");
        } else {
            $width *= $scale;
            $height *= $scale;
	        open (PNG, '|-', "rsvg --format=png -w ${width} -h ${height} $filename /dev/stdout") or error("Couldn't convert to PNG.");
        }
	    print <PNG>;
	    close(PNG) or error("Failed to close pipe.");
    }
    else {
	    print $output;
    }
}
=back
=cut

1;
