
sub do_env_classinterface {
    local ($_) = @_;
    local ($list_type) = 'DL';

    s/$next_pair_rx//;     # Ditch the label specifier
    s/$next_pair_rx//;     # Ditto the length declarations ...
                           # but we may want to switch to enumerated style
                           # if they include a \usecounter.

    $list_type = 'OL' if $1 =~ /\\usecounter/;

    &list_helper($_, $list_type);
}

sub do_cmd_clsnm {
    local($_) = @_;
    local($rest) = $_;
    $rest =~ s/$next_pair_pr_rx//o;
    local($match) = $&;
    join('',"<tt>",$&,"<\/tt>",$rest);
}

sub do_cmd_clsnmref {
    local($_) = @_;
    local($rest) = $_;
    $rest =~ s/$next_pair_pr_rx//o;
    local($match) = $&;
    &process_ref($cross_ref_mark,$cross_ref_mark,$match);
}

sub do_cmd_vrbl {
    local($_) = @_;
    local($rest) = $_;
    $rest =~ s/$next_pair_pr_rx//o;
    local($match) = $&;
    join('',"<i>",$&,"<\/i>",$rest);
}

sub do_cmd_srccd {
    local($_) = @_;
    local($rest) = $_;
    $rest =~ s/$next_pair_pr_rx//o;
    local($match) = $&;
    join('',"<tt>",$&,"<\/tt>",$rest);
}

sub do_cmd_exenm {
    local($_) = @_;
    local($rest) = $_;
    $rest =~ s/$next_pair_pr_rx//o;
    local($match) = $&;
    join('',"<tt>",$&,"<\/tt>",$rest);
}

sub do_cmd_type {
    local($_) = @_;
    local($rest) = $_;
    $rest =~ s/$next_pair_pr_rx//o;
    local($match) = $&;
    join('',"<tt>",$&,"<\/tt>",$rest);
}

sub do_cmd_keywd {
    local($_) = @_;
    local($rest) = $_;
    $rest =~ s/$next_pair_pr_rx//o;
    local($match) = $&;
    join('',"<tt>",$&,"<\/tt>",$rest);
}

sub do_cmd_filnm {
    local($_) = @_;
    local($rest) = $_;
    $rest =~ s/$next_pair_pr_rx//o;
    local($match) = $&;
    join('',"<tt>",$&,"<\/tt>",$rest);
}

sub do_cmd_clssection {
    local($_) = @_;
    s/$next_pair_pr_rx//o;
    local($label) = $2;
    local($header) = "<h3>The <tt>$label</tt> Class</h3>";
    join('', $header, &anchor_label($label,$CURRENT_FILE,$_));
}

sub do_cmd_clssubsection {
    local($_) = @_;
    local($rest) = $_;
    $rest =~ s/$next_pair_pr_rx//o;
    local($match) = $&;
    join('',"<h4>",$&,"<\/h4>",$rest);
}

1;                              # This must be the last line
