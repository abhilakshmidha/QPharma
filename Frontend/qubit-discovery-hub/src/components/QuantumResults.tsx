import React, { useState } from "react";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Progress } from "@/components/ui/progress";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import {
  Atom,
  Download,
  TrendingUp,
  Zap,
  Target,
  ChevronRight,
  AlertCircle,
  CheckCircle,
  Clock,
  FileText,
  Eye,
  BarChart3,
} from "lucide-react";
import { cn } from "@/lib/utils";
import { toast } from "sonner";

interface QuantumResultsProps {
  results: {
    session_id: string;
    status: string;
    message: string;
    results: {
      classical?: Record<string, any>;
      quantum?: Record<string, any>;
      processing_times: {
        classical?: Record<string, number>;
        quantum?: Record<string, number>;
      };
    };
  };
}

export const QuantumResults: React.FC<QuantumResultsProps> = ({ results }) => {
  const [summaryData, setSummaryData] = useState(null);
  const [loadingSummary, setLoadingSummary] = useState(false);

  const hasClassical =
    results.results.classical &&
    Object.keys(results.results.classical).length > 0;
  const hasQuantum =
    results.results.quantum && Object.keys(results.results.quantum).length > 0;

  const fetchSummary = async () => {
    setLoadingSummary(true);
    try {
      const response = await fetch(`/api/summary/${results.session_id}`);
      if (response.ok) {
        const data = await response.json();
        setSummaryData(data);
        toast.success("Summary generated successfully!");
      } else {
        throw new Error("Failed to fetch summary");
      }
    } catch (error) {
      toast.error("Failed to generate summary");
    } finally {
      setLoadingSummary(false);
    }
  };

  const downloadResults = async () => {
    try {
      const response = await fetch(`/api/results/${results.session_id}`);
      if (response.ok) {
        const data = await response.json();
        const blob = new Blob([JSON.stringify(data, null, 2)], {
          type: "application/json",
        });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = `results_${results.session_id}.json`;
        a.click();
        URL.revokeObjectURL(url);
        toast.success("Results downloaded successfully!");
      }
    } catch (error) {
      toast.error("Failed to download results");
    }
  };

  const openVisualization = (htmlContent: string, title: string) => {
    try {
      const decodedHtml = atob(htmlContent);
      const newWindow = window.open("", "_blank");
      if (newWindow) {
        newWindow.document.write(decodedHtml);
        newWindow.document.title = title;
        newWindow.document.close();
      }
    } catch (error) {
      toast.error("Failed to open visualization");
    }
  };

  const renderFileResults = (
    fileResults: any,
    type: "classical" | "quantum"
  ) => {
    return Object.entries(fileResults).map(
      ([filename, result]: [string, any]) => (
        <Card
          key={`${type}-${filename}`}
          className="bg-gradient-card border-border/50"
        >
          <CardHeader className="pb-3">
            <div className="flex items-center justify-between">
              <CardTitle className="text-lg flex items-center space-x-2">
                <FileText className="w-5 h-5" />
                <span>{filename}</span>
              </CardTitle>
              <Badge variant={result.success ? "default" : "destructive"}>
                {result.success ? "Success" : "Failed"}
              </Badge>
            </div>
          </CardHeader>

          <CardContent className="space-y-4">
            {result.success ? (
              <>
                {/* Processing Time */}
                <div className="flex items-center justify-between">
                  <div className="flex items-center space-x-2">
                    <Clock className="w-4 h-4 text-muted-foreground" />
                    <span className="text-sm text-muted-foreground">
                      Processing Time
                    </span>
                  </div>
                  <span className="font-mono font-medium">
                    {result.processing_time?.toFixed(4) || "N/A"}s
                  </span>
                </div>

                {/* Molecule Info */}
                {result.molecule_info && (
                  <div className="space-y-2">
                    <h4 className="font-semibold text-sm">
                      Molecular Properties
                    </h4>
                    <div className="grid grid-cols-2 gap-2 text-sm">
                      <div className="flex justify-between">
                        <span className="text-muted-foreground">Atoms:</span>
                        <span>{result.molecule_info.num_atoms}</span>
                      </div>
                      <div className="flex justify-between">
                        <span className="text-muted-foreground">Bonds:</span>
                        <span>{result.molecule_info.num_bonds}</span>
                      </div>
                      <div className="flex justify-between col-span-2">
                        <span className="text-muted-foreground">
                          Molecular Weight:
                        </span>
                        <span>
                          {result.molecule_info.molecular_weight?.toFixed(3) ||
                            "N/A"}{" "}
                          g/mol
                        </span>
                      </div>
                    </div>
                  </div>
                )}

                {/* Quantum Output */}
                {type === "quantum" && result.quantum_output && (
                  <div className="space-y-2">
                    <h4 className="font-semibold text-sm">
                      Quantum Simulation Output
                    </h4>
                    <div className="bg-muted/30 p-3 rounded text-xs font-mono">
                      {JSON.stringify(result.quantum_output, null, 2)}
                    </div>
                  </div>
                )}

                {/* Visualization Button */}
                {result.visualization && (
                  <Button
                    variant="outline"
                    size="sm"
                    onClick={() =>
                      openVisualization(
                        result.visualization,
                        `${type} - ${filename}`
                      )
                    }
                    className="w-full"
                  >
                    <Eye className="w-4 h-4 mr-2" />
                    View 3D Visualization
                  </Button>
                )}
              </>
            ) : (
              <div className="flex items-center space-x-2 text-destructive">
                <AlertCircle className="w-4 h-4" />
                <span className="text-sm">
                  {result.error || "Processing failed"}
                </span>
              </div>
            )}
          </CardContent>
        </Card>
      )
    );
  };

  return (
    <div className="w-full max-w-6xl mx-auto space-y-6">
      {/* Header */}
      <Card className="bg-gradient-card border-border/50">
        <CardHeader className="pb-4">
          <div className="flex items-center justify-between">
            <div className="flex items-center space-x-3">
              <div className="p-2 rounded-lg bg-gradient-quantum">
                <Atom className="w-6 h-6 text-primary-foreground" />
              </div>
              <div>
                <CardTitle className="text-xl">Analysis Results</CardTitle>
                <p className="text-muted-foreground">
                  Session: {results.session_id}
                </p>
              </div>
            </div>
            <Badge variant="default" className="animate-quantum-pulse">
              {results.status}
            </Badge>
          </div>
        </CardHeader>
        <CardContent>
          <p className="text-muted-foreground">{results.message}</p>
        </CardContent>
      </Card>

      {/* Processing Times Summary */}
      {(hasClassical || hasQuantum) && (
        <Card className="bg-gradient-card border-border/50">
          <CardHeader>
            <CardTitle className="text-lg flex items-center space-x-2">
              <TrendingUp className="w-5 h-5" />
              <span>Processing Performance</span>
            </CardTitle>
          </CardHeader>
          <CardContent>
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
              {hasClassical && (
                <div className="space-y-2">
                  <div className="flex items-center space-x-2">
                    <Zap className="w-4 h-4 text-blue-500" />
                    <span className="font-medium">Classical Processing</span>
                  </div>
                  <div className="space-y-1">
                    {Object.entries(
                      results.results.processing_times.classical || {}
                    ).map(([file, time]) => (
                      <div key={file} className="flex justify-between text-sm">
                        <span className="text-muted-foreground truncate">
                          {file}
                        </span>
                        <span className="font-mono">
                          {(time as number).toFixed(4)}s
                        </span>
                      </div>
                    ))}
                  </div>
                </div>
              )}

              {hasQuantum && (
                <div className="space-y-2">
                  <div className="flex items-center space-x-2">
                    <Target className="w-4 h-4 text-purple-500" />
                    <span className="font-medium">Quantum Processing</span>
                  </div>
                  <div className="space-y-1">
                    {Object.entries(
                      results.results.processing_times.quantum || {}
                    ).map(([file, time]) => (
                      <div key={file} className="flex justify-between text-sm">
                        <span className="text-muted-foreground truncate">
                          {file}
                        </span>
                        <span className="font-mono">
                          {(time as number).toFixed(4)}s
                        </span>
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </div>
          </CardContent>
        </Card>
      )}

      {/* Results Tabs */}
      <Tabs
        defaultValue={hasClassical ? "classical" : "quantum"}
        className="w-full"
      >
        <TabsList className="grid w-full grid-cols-2">
          {hasClassical && (
            <TabsTrigger
              value="classical"
              className="flex items-center space-x-2"
            >
              <Zap className="w-4 h-4" />
              <span>Classical Results</span>
            </TabsTrigger>
          )}
          {hasQuantum && (
            <TabsTrigger
              value="quantum"
              className="flex items-center space-x-2"
            >
              <Target className="w-4 h-4" />
              <span>Quantum Results</span>
            </TabsTrigger>
          )}
        </TabsList>

        {hasClassical && (
          <TabsContent value="classical" className="space-y-4">
            <div className="grid gap-4">
              {renderFileResults(results.results.classical, "classical")}
            </div>
          </TabsContent>
        )}

        {hasQuantum && (
          <TabsContent value="quantum" className="space-y-4">
            <div className="grid gap-4">
              {renderFileResults(results.results.quantum, "quantum")}
            </div>
          </TabsContent>
        )}
      </Tabs>

      {/* Summary Section */}
      {summaryData && (
        <Card className="bg-gradient-card border-border/50">
          <CardHeader>
            <CardTitle className="text-lg flex items-center space-x-2">
              <BarChart3 className="w-5 h-5" />
              <span>Performance Summary</span>
            </CardTitle>
          </CardHeader>
          <CardContent>
            {summaryData.summary_plots?.summary_plot && (
              <div className="mb-4">
                <img
                  src={`data:image/png;base64,${summaryData.summary_plots.summary_plot}`}
                  alt="Performance Summary"
                  className="w-full max-w-4xl mx-auto rounded-lg border"
                />
              </div>
            )}

            {summaryData.summary_plots?.statistics && (
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-sm">
                {Object.entries(summaryData.summary_plots.statistics).map(
                  ([method, stats]: [string, any]) => (
                    <div key={method} className="space-y-2">
                      <h4 className="font-semibold capitalize">
                        {method} Statistics
                      </h4>
                      <div className="space-y-1 text-muted-foreground">
                        <div className="flex justify-between">
                          <span>Total Time:</span>
                          <span>{stats.total_time?.toFixed(4)}s</span>
                        </div>
                        <div className="flex justify-between">
                          <span>Average Time:</span>
                          <span>{stats.average_time?.toFixed(4)}s</span>
                        </div>
                        <div className="flex justify-between">
                          <span>Min Time:</span>
                          <span>{stats.min_time?.toFixed(4)}s</span>
                        </div>
                        <div className="flex justify-between">
                          <span>Max Time:</span>
                          <span>{stats.max_time?.toFixed(4)}s</span>
                        </div>
                      </div>
                    </div>
                  )
                )}
              </div>
            )}
          </CardContent>
        </Card>
      )}

      {/* Action Buttons */}
      <div className="flex flex-wrap gap-4 justify-center">
        <Button variant="outline" size="lg" onClick={downloadResults}>
          <Download className="w-4 h-4 mr-2" />
          Download Results
        </Button>
        <Button
          variant="quantum"
          size="lg"
          onClick={fetchSummary}
          disabled={loadingSummary}
        >
          <BarChart3 className="w-4 h-4 mr-2" />
          {loadingSummary ? "Generating..." : "Generate Summary"}
        </Button>
      </div>
    </div>
  );
};
