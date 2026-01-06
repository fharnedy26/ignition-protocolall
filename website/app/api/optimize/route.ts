import { NextRequest, NextResponse } from "next/server";
import { exec } from "child_process";
import { promisify } from "util";
import { writeFile, unlink, mkdir, readFile, readdir } from "fs/promises";
import { join } from "path";
import { existsSync } from "fs";

const execAsync = promisify(exec);

// DEMO MODE: This API route is currently in demo mode.
// For production deployment, implement full engine integration by:
// 1. Parsing the Top-N.csv output file from engine.py
// 2. Reading the RUN.json metadata file
// 3. Returning actual optimization results
// 
// Note: The engine.py integration requires Python environment and dependencies.
// For expo/demo purposes, this route returns simulated results.

const DEMO_MODE = true; // Set to false when implementing full integration

export async function POST(request: NextRequest) {
  try {
    const body = await request.json();
    const { components, fuel_type = "gasoline", K = 5, top = 10 } = body;

    if (!components || !Array.isArray(components) || components.length === 0) {
      return NextResponse.json(
        { error: "Components array is required" },
        { status: 400 }
      );
    }

    // Validate fuel_type
    const validFuelTypes = ["gasoline", "diesel", "jet"];
    if (!validFuelTypes.includes(fuel_type.toLowerCase())) {
      return NextResponse.json(
        { error: `Invalid fuel_type. Must be one of: ${validFuelTypes.join(", ")}` },
        { status: 400 }
      );
    }

    // Create temporary CSV file
    const tempDir = join(process.cwd(), "tmp");
    if (!existsSync(tempDir)) {
      await mkdir(tempDir, { recursive: true });
    }

    const timestamp = Date.now();
    const csvPath = join(tempDir, `components_${timestamp}.csv`);
    
    // Convert components to CSV
    const headers = Object.keys(components[0]);
    const csvContent = [
      headers.join(","),
      ...components.map((comp: any) =>
        headers.map((h) => {
          const val = comp[h];
          // Escape commas and quotes in CSV
          if (typeof val === "string" && (val.includes(",") || val.includes('"'))) {
            return `"${val.replace(/"/g, '""')}"`;
          }
          return val ?? "";
        }).join(",")
      ),
    ].join("\n");

    await writeFile(csvPath, csvContent);

    try {
      if (DEMO_MODE) {
        // DEMO MODE: Return simulated results
        // This allows the demo to work without requiring Python engine setup
        return NextResponse.json({
          success: true,
          demo: true,
          message: "Optimization completed (DEMO MODE - Simulated Results)",
          warning: "This is a demonstration. For production use, implement full engine integration.",
          score: 85.5 + Math.random() * 10,
          properties: {
            RON: 90 + Math.random() * 10,
            MON: 85 + Math.random() * 10,
            LHV: 40 + Math.random() * 5,
            PMI: 10 + Math.random() * 10,
            cetane: fuel_type === "diesel" ? 45 + Math.random() * 10 : 0,
          },
          components: components.slice(0, parseInt(String(K))).map((comp: any, idx: number) => ({
            ...comp,
            vol_frac: (1 / parseInt(String(K))) + (Math.random() - 0.5) * 0.1,
          })),
        });
      }

      // PRODUCTION MODE: Call actual engine.py
      const enginePath = join(process.cwd(), "..", "engine.py");
      const outputDir = join(tempDir, `output_${timestamp}`);

      // Ensure output directory exists
      if (!existsSync(outputDir)) {
        await mkdir(outputDir, { recursive: true });
      }

      const command = `python "${enginePath}" --components "${csvPath}" --fuel-type ${fuel_type} --K ${K} --top ${top} --seed 42 --out "${outputDir}"`;

      try {
        const { stdout, stderr } = await execAsync(command, {
          cwd: join(process.cwd(), ".."),
          timeout: 30000, // 30 second timeout
        });

        // Parse results from Top-N.csv
        const topNFile = join(outputDir, `Top-${top}.csv`);
        if (existsSync(topNFile)) {
          const csvContent = await readFile(topNFile, "utf-8");
          const lines = csvContent.split("\n").filter((line) => line.trim());
          const headers = lines[0].split(",");
          
          if (lines.length > 1) {
            const bestResult = lines[1].split(",");
            const result: any = {};
            headers.forEach((header, idx) => {
              result[header.trim()] = bestResult[idx]?.trim();
            });

            // Parse component fractions from the components column
            const componentsStr = result.components || "";
            const componentList = componentsStr.split(";").map((comp: string) => {
              const [id, frac] = comp.split(":");
              return {
                id: id?.trim(),
                vol_frac: parseFloat(frac?.trim() || "0"),
              };
            });

            return NextResponse.json({
              success: true,
              demo: false,
              message: "Optimization completed",
              score: parseFloat(result.score || "0"),
              properties: {
                RON: parseFloat(result.RON || "0"),
                MON: parseFloat(result.MON || "0"),
                LHV: parseFloat(result.LHV_mass || "0"),
                PMI: parseFloat(result.PMI || "0"),
                cetane: parseFloat(result.cetane || "0"),
              },
              components: componentList,
              raw: result, // Include full result for debugging
            });
          }
        }

        // Fallback if parsing fails
        return NextResponse.json({
          success: true,
          demo: false,
          message: "Optimization completed but result parsing failed",
          warning: "Check server logs for details",
          stdout: stdout.substring(0, 500), // First 500 chars for debugging
        });
      } catch (execError: any) {
        // Engine execution failed - return error
        console.error("Engine execution error:", execError);
        return NextResponse.json(
          {
            error: "Engine execution failed",
            message: execError.message || "Unknown error",
            details: process.env.NODE_ENV === "development" ? execError.stderr : undefined,
            demo: false,
          },
          { status: 500 }
        );
      }
    } finally {
      // Cleanup temporary files
      try {
        if (existsSync(csvPath)) {
          await unlink(csvPath);
        }
      } catch (e) {
        // Ignore cleanup errors
      }
    }
  } catch (error: any) {
    console.error("Optimization error:", error);
    return NextResponse.json(
      {
        error: error.message || "Optimization failed",
        details: process.env.NODE_ENV === "development" ? error.stack : undefined,
      },
      { status: 500 }
    );
  }
}

