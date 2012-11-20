# Converts the output of git diff -C <commit1> <commit2> or git show -C <commit>
# to html
#
# Author: David Roazen
# Modified: 7/11/2011
#

use strict;
use Getopt::Long qw(GetOptions);

use constant TRUE  => 1;
use constant FALSE => 0;

my $NORMAL_DIFF_DELIMITER = "diff --git";
my $MERGE_DIFF_DELIMITER = "diff --cc";
my $PRE_OPTIONS = "style=\"font-size: 115%; white-space: pre-wrap;\"";
my $div_background_color = "eaf2f5";

my @print_buffer;
my @table_of_contents_links;
my @table_of_contents_files;
my $table_of_contents_header = "Changes";
my $commit_header_style = "anchor";
my $section_number = 0;
my $commit_id;
my $old_revision;
my $new_revision;
my $diff_only = FALSE;

GetOptions("toc_header=s" => \$table_of_contents_header,
           "commit_header_style=s" => \$commit_header_style,
           "old_revision=s" => \$old_revision,
           "new_revision=s" => \$new_revision,
           "div_color=s" => \$div_background_color
           );            
undef @ARGV;

convert();
exit(0);

sub convert {    
    my $first_line = <>;
    my $next_diff_start = "";
    
    unless ( is_diff_start($first_line) or is_commit_header_start($first_line) ) {
        echo_all_input_verbatim($first_line);
        exit(1);
    }

    print "<div width=\"80%\" style=\"background-color:#${div_background_color};\">\n";
    
    if ( is_commit_header_start($first_line) ) {
        $next_diff_start = convert_commit_header($first_line);            
    }
    else {
        $next_diff_start = $first_line;
        $diff_only = TRUE;
    }
    
    while ( $next_diff_start ) {
        $next_diff_start = convert_diff($next_diff_start);
    }
    
    print_table_of_contents();
    flush_print_buffer();
    
    print "</div><br />\n";
}

sub convert_commit_header {
    my $commit_line = shift; chomp $commit_line;
        
    ($commit_id) = $commit_line =~ /commit (\w+)/;
    
    escape($commit_line, $commit_id);
    
    if ( $commit_header_style eq "link" ) {
        printf "<a href=\"#%s\"><b>%s</b></a><br />\n", $commit_id, $commit_line;
    }
    elsif ( $commit_header_style eq "anchor" ) {
        printf "<a name=\"%s\"><b>%s</b></a><br />\n", $commit_id, $commit_line;
    }
    else {
        printf "<b>%s</b><br />\n", $commit_line;
    }

    my $header_line;
    while ( $header_line = <> ) {
        last if ( $header_line !~ /^\w+\:.*/ );
        
        escape($header_line);
        $header_line =~ s/^(\w+\:)(.*)/<b>$1<\/b>$2/;
        
        print "${header_line}<br />\n";
    }
    
    return convert_commit_header_log_message($header_line);
}

sub convert_commit_header_log_message {
    my $message_line = shift;
    
    return $message_line unless $message_line;
    
    print "<pre ${PRE_OPTIONS}>"; 
    do {
        if ( is_diff_start($message_line) ) {
            print "</pre>\n";
            return $message_line;
        }
        
        escape($message_line);
        print $message_line;
    } while ( $message_line = <> );
    
    print "</pre>\n";   
    return $message_line;
}

sub convert_diff {
    my $last_diff_summary_line = convert_diff_summary(shift);
    
    if ( is_diff_start($last_diff_summary_line) ) {
        return $last_diff_summary_line;
    }
    
    return convert_diff_text();
}

sub convert_diff_summary {
    my $diff_start = shift; 
    chomp $diff_start;
    
    my @diff_header_lines;
    push @diff_header_lines, $diff_start;

    my $action = "Modified";
    my $line;
    
    while ( $line = <> ) {
        chomp $line;
                
        if ( $line =~ /^deleted file mode/ ) {
            $action = "Deleted";
        }
        elsif ( $line =~ /^new file mode/ ) {
            $action = "Created";
        }
        elsif ( $line =~ /^rename from/ ) {
            $action = "Renamed";
        }
        elsif ( $line =~ /^copy from/ ) {
            $action = "Copied";
        }
        
        unless ( is_diff_start($line) ) {
            push @diff_header_lines, $line;
        }
        
        last if $line =~ /^index/ or is_diff_start($line);
    }

    my ($file_name, $new_file_name);
    if ( $diff_start =~ /^${NORMAL_DIFF_DELIMITER}/ ) {
        ($file_name, $new_file_name) = 
            $diff_start =~ /^${NORMAL_DIFF_DELIMITER} a\/(\S+) b\/(\S+)/;
        escape($new_file_name);
    }
    elsif ( $diff_start =~ /^${MERGE_DIFF_DELIMITER}/ ) {
        ($file_name) = $diff_start =~ /^${MERGE_DIFF_DELIMITER} (\S+)/;
    }
        
    escape($file_name);
    
    my $section_id = get_new_section_id();
    
    my $section_heading = "$action file: $file_name";
    if ( $action eq "Renamed" or $action eq "Copied" ) {
        $section_heading .= " to $new_file_name";
        push @table_of_contents_files, "$file_name $new_file_name";
    }
    else {
        push @table_of_contents_files, $file_name;
    }
    
    push @table_of_contents_links, "<a href=\"#${section_id}\">${section_heading}</a><br />";
    
    print_to_buffer("<a name=\"${section_id}\"><b>${section_heading}</b></a>\n");   
    print_to_buffer("<hr />");
        
    print_to_buffer("<pre ${PRE_OPTIONS}>");
    for my $diff_header_line ( @diff_header_lines ) {
        escape($diff_header_line);
        print_to_buffer("$diff_header_line\n");
    }
    print_to_buffer("</pre>");
    
    return $line;
}

sub convert_diff_text {
    my $line;

    print_to_buffer("<pre ${PRE_OPTIONS}>");

    while ( $line = <> ) {        
        last if is_diff_start($line);

        if ( $line =~ /^@@.+?@@/ ) {
            chomp $line;
            escape($line);
            print_to_buffer("</pre>\n");
            print_to_buffer("<i>${line}</i>\n");
            print_to_buffer("<pre ${PRE_OPTIONS}>");
        }
        else {
            escape($line);
            print_to_buffer($line);
        }
    }
    print_to_buffer("</pre>\n");

    return $line;
}

sub escape {
    for my $i ( 0..$#_ ) {
        $_[$i] =~ s/&/&amp;/g;
        $_[$i] =~ s/</&lt;/g;
        $_[$i] =~ s/>/&gt;/g;        
    }
}

sub is_commit_header_start {
    my $candidate = shift;
    
    if ( $candidate =~ /^commit \w+/ ) {
        return TRUE;
    }
    
    return FALSE;
}

sub is_diff_start {
    my $candidate = shift;
    
    if ( $candidate =~ /^${NORMAL_DIFF_DELIMITER}/ or 
         $candidate =~ /^${MERGE_DIFF_DELIMITER}/ ) {
        return TRUE;
    }
    
    return FALSE;
}

sub get_new_section_id {
    my $section_id;
    
    if ( $commit_id ) {
        $section_id = "${commit_id}_${section_number}";
    }
    else {
        $section_id = "diff_${section_number}";
    }
    
    $section_number++;
    return $section_id;
}

sub print_to_buffer {
    my $line = shift;
    
    push @print_buffer, $line;
}

sub flush_print_buffer {
    for my $line ( @print_buffer ) {
        print $line;
    }
}

sub print_table_of_contents {
    return unless @table_of_contents_links and @table_of_contents_files;
    
    if ( $diff_only ) {
        print "<a name=\"diff\"><b>${table_of_contents_header}:</b></a><hr />";
    }
    else {
        print "<b>${table_of_contents_header}:</b><hr />";
    }
    
    if ( $old_revision and $new_revision ) {
        print "<table border=\"0\" cellpadding=\"3\" cellspacing=\"0\">\n";
        
        for my $i ( 0..$#table_of_contents_links ) {
            my $toc_entry = $table_of_contents_links[$i];
            my $toc_file = $table_of_contents_files[$i];
            
            my $file_numstat = `git diff-tree -C --numstat $old_revision $new_revision -- $toc_file`;
            my ($file_additions, $file_deletions) = 
                $file_numstat =~ /^(\S+)\s+(\S+)\s+/;
            
            if ( $file_additions =~ /^\-$/ ) {
                $file_additions = "(binary)";
            }
            else {
                $file_additions .= "+";
            }
            
            if ( $file_deletions =~ /^\-$/ ) {
                $file_deletions = "(binary)";
            }
            else {
                $file_deletions .= "-";
            }            
            
            printf "<tr><td><b>%s</b></td><td><b>%s</b></td><td><b>%s</b></td></tr>\n",
                   $file_additions, $file_deletions, $toc_entry; 
        }
        
        print "</table>\n";
    }
    else {
        for my $toc_entry ( @table_of_contents_links ) {
            print "${toc_entry}\n";
        }
    }
    print "<br /><br />";
}

sub echo_all_input_verbatim {
    my $line = shift;
    
    print "<pre ${PRE_OPTIONS}>\n";
    do {
        escape($line);
        print $line;
    } while ( $line = <> );
    print "</pre>\n";
}