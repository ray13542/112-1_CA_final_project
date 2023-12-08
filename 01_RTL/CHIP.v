//----------------------------- DO NOT MODIFY THE I/O INTERFACE!! ------------------------------//
module CHIP #(                                                                                  //
    parameter BIT_W = 32                                                                        //
)(                                                                                              //
    // clock                                                                                    //
        input               i_clk,                                                              //
        input               i_rst_n,                                                            //
    // instruction memory                                                                       //
        input  [BIT_W-1:0]  i_IMEM_data,                                                        //
        output [BIT_W-1:0]  o_IMEM_addr,                                                        //
        output              o_IMEM_cen,                                                         //
    // data memory                                                                              //
        input               i_DMEM_stall,                                                       //
        input  [BIT_W-1:0]  i_DMEM_rdata,                                                       //
        output              o_DMEM_cen,                                                         //
        output              o_DMEM_wen,                                                         //
        output [BIT_W-1:0]  o_DMEM_addr,                                                        //
        output [BIT_W-1:0]  o_DMEM_wdata,                                                       //
    // finnish procedure                                                                        //
        output              o_finish,                                                           //
    // cache                                                                                    //
        input               i_cache_finish,                                                     //
        output              o_proc_finish                                                       //
);                                                                                              //
//----------------------------- DO NOT MODIFY THE I/O INTERFACE!! ------------------------------//

// ------------------------------------------------------------------------------------------------------------------------------------------------------
// Parameters
// ------------------------------------------------------------------------------------------------------------------------------------------------------

    // TODO: any declaration
    // func7 opcode
    parameter R_type = 7'b0110011; // add, sub, and, xor, mul
    parameter I_type = 7'b0010011; // addi, slli, slti, srai
    parameter auipc_type  = 7'b0010111; // auipc
    parameter sw_type = 7'b0100011; // sw
    parameter lw_type = 7'b0000011; // lw
    parameter B_type = 7'b1100011; // beq, bge, blt, bne
    parameter jal_type = 7'b1101111; // jal
    parameter jalr_type = 7'b1100111; // jalr
    parameter ecall_type = 7'b1110011; // ecall

// ------------------------------------------------------------------------------------------------------------------------------------------------------
// Wires and Registers
// ------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // TODO: any declaration
    reg [BIT_W-1:0] PC, next_PC;
    reg [BIT_W-1:0] PCadd4, PCbranch;

    wire mem_cen, mem_wen;
    wire [BIT_W-1:0] mem_addr, mem_wdata, mem_rdata;
    wire mem_stall;
    wire [6:0] opcode;
    wire [2:0] funct3;
    wire [4:0] rs1, rs2, rd;

    wire [BIT_W-1:0] reg_rdata_1, reg_rdata_2, alu_in1, alu_in2, alu_result;
    //control signal
    wire zero;
    wire goBranch, Branch, MemRead, MemtoReg, MemWrite, ALUsrc, RegWrite;
    wire [1:0] ALUop;
    wire [3:0] ALUctrl;
        


// ------------------------------------------------------------------------------------------------------------------------------------------------------
// Continuous Assignment
// ------------------------------------------------------------------------------------------------------------------------------------------------------
    // TODO: any wire assignment
    assign PC = o_IMEM_addr;
    assign opcode = i_IMEM_data[6:0];
    assign funct3 = i_IMEM_data[14:12];
    assign rs2 = i_IMEM_data[24:20];
    assign rs1 = i_IMEM_data[19:15];
    assign rd = i_IMEM_data[11:7];
    assign goBranch = Branch & zero;
    assign alu_in1 = reg_rdata_1;
    assign alu_in2 = (ALUsrc)? imm : reg_rdata_2;
// ------------------------------------------------------------------------------------------------------------------------------------------------------
// Submoddules
// ------------------------------------------------------------------------------------------------------------------------------------------------------
    ALU alu(.in1(alu_in1), .in2(alu_in2), .ALUctrl(ALUctrl), .result(alu_result), zero.(zero));
    MULDIV_unit mul(.result(aluresult), .o_done(), .i_clk(i_clk), .i_valid(), .i_A(alu_in1), .i_B(alu_in2), .ALUctrl(ALUctrl));
    ALUcontrol alucontrol(.ALUop(ALUop), .ALUctrl(ALUctrl), .funct7(), .funct3(funct3));
    // TODO: Reg_file wire connection
    Reg_file reg0(               
        .i_clk  (i_clk),             
        .i_rst_n(i_rst_n),         
        .wen    (RegWrite),          
        .rs1    (rs1),                
        .rs2    (rs2),                
        .rd     (rd),                 
        .wdata  (),             
        .rdata1 (reg_rdata_1),           
        .rdata2 (reg_rdata_2)
    );

// ------------------------------------------------------------------------------------------------------------------------------------------------------
// Always Blocks
// ------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Todo: any combinational/sequential circuit
    
    // Control 
    // reg Branch, MemRead, MemtoReg, MemWrite, ALUsrc, RegWrite;
    reg jal,jalr;
    always @(*) begin
        ALUsrc = 0;
        RegWrite = 0;
        Branch = 0;
        MemtoReg = 0;
        MemWrite = 0;
        MemRead = 0;
        jal = 0;
        jalr = 0;
        // 
        case(opcode)
            R_type: begin
                ALUsrc = 0;
                RegWrite = 1;
                ALUop = 2'b10;
                MemWrite = 0;
                MemRead = 0;
                MemtoReg = 0;
            end
            I_type: begin // addi, slli, slti, srai
                ALUsrc = 1;
                ALUop = 2'b11;
                RegWrite = 1;
                MemWrite = 0;
                MemRead = 0;
                MemtoReg = 0;
            end
            auipc_type: begin
                ALUsrc = 1;
                ALUop = 2'b0;
                RegWrite = 1; // rd stores pc+imm
                MemWrite = 0;
                MemRead = 0;
                MemtoReg = 0; //
            end
            sw_type: begin
                ALUsrc = 1;
                ALUop = 2'b0;  
                RegWrite = 0;
                MemWrite = 1;
                MemRead = 0;
                MemtoReg = 0;
            end
            lw_type: begin
                ALUsrc = 1;
                ALUop = 2'b0;
                RegWrite = 1;
                MemWrite = 0;
                MemRead = 1;
                MemtoReg = 1;
            end
            B_type: begin
                ALUsrc = 0;
                ALUop = 2'b01;
                RegWrite = 0;
                MemWrite = 0;
                MemRead = 0;
                MemtoReg = 0;
                Branch = 1; // branch
            end
            jal_type: begin
                ALUsrc = 1;
                ALUop = 2'b0;
                RegWrite = 1; // write pc+4 to rd
                MemWrite = 0;
                MemRead = 0;
                MemtoReg = 0; 
            end
            jalr_type: begin
                ALUsrc = 1;
                ALUop = 2'b0;
                RegWrite = 1; // write pc+4 to rd
                MemWrite = 0;
                MemRead = 0;
                MemtoReg = 0;
            end
        endcase
    end
    
    // ImmGen part
    reg [31:0] imm;
    always@(*) begin
        imm = 0;
        case(opcode)
            I_type: begin // addi, slli, slti, srai
                if(funct3 == 3'b001 || funct3 == 3'b101) // slli, srai are different 
                    imm = {{27{i_IMEM_data[24]}}, i_IMEM_data[24:20]};
                else // addi, slti
                    imm = {{20{i_IMEM_data[31]}}, i_IMEM_data[31:20]}; // 12 bits imm, sign extension
            end
            auipc_type: 
                imm = {i_IMEM_data[31:12], 12'b0};
            sw_type: 
                imm = {{20{i_IMEM_data[31]}}, i_IMEM_data[31:25], i_IMEM_data[11:7]};
            lw_type: // lw is I-type
                imm = {{20{i_IMEM_data[31]}}, i_IMEM_data[31:20]}; 
            beq_type: 
                imm = {{20{i_IMEM_data[31]}}, i_IMEM_data[7], i_IMEM_data[30:25], i_IMEM_data[11:8], 1'b0}; 
            jal_type:
                imm = {{12{i_IMEM_data[31]}}, i_IMEM_data[19:12], i_IMEM_data[20], i_IMEM_data[30:21], 1'b0};
            jalr_type: // jalr is I-type
                imm = {{20{i_IMEM_data[31]}}, i_IMEM_data[31:20]};
        endcase
    end
    // OK -> still need to do : When ALUsrc = 1, imm passes to ALU ; When ALUsrc = 0, Read_data_2 passes to ALU
                
    
    always @(posedge i_clk) begin
        PCadd4 = PC + 4;
        PCbranch = (imm << 1) + PC;
        next_PC = ()? PCbranch: PCadd4;
    end

    always @(posedge i_clk or negedge i_rst_n) begin
        if (!i_rst_n) begin
            PC <= 32'h00010000; // Do not modify this value!!!
        end
        else begin
            PC <= next_PC;
        end
    end
endmodule

module Reg_file(i_clk, i_rst_n, wen, rs1, rs2, rd, wdata, rdata1, rdata2);
   
    parameter BITS = 32;
    parameter word_depth = 32;
    parameter addr_width = 5; // 2^addr_width >= word_depth
    
    input i_clk, i_rst_n, wen; // wen: 0:read | 1:write
    input [BITS-1:0] wdata;
    input [addr_width-1:0] rs1, rs2, rd;

    output [BITS-1:0] rdata1, rdata2;

    reg [BITS-1:0] mem [0:word_depth-1];
    reg [BITS-1:0] mem_nxt [0:word_depth-1];

    integer i;

    assign rdata1 = mem[rs1];
    assign rdata2 = mem[rs2];

    always @(*) begin
        for (i=0; i<word_depth; i=i+1)
            mem_nxt[i] = (wen && (rd == i)) ? wdata : mem[i];
    end

    always @(posedge i_clk or negedge i_rst_n) begin
        if (!i_rst_n) begin
            mem[0] <= 0;
            for (i=1; i<word_depth; i=i+1) begin
                case(i)
                    32'd2: mem[i] <= 32'hbffffff0;
                    32'd3: mem[i] <= 32'h10008000;
                    default: mem[i] <= 32'h0;
                endcase
            end
        end
        else begin
            mem[0] <= 0;
            for (i=1; i<word_depth; i=i+1)
                mem[i] <= mem_nxt[i];
        end       
    end
endmodule
module ALUcontrol(ALUop, ALUctrl, funct7, funct3);
    input [1:0] ALUop;
    input [2:0] funct3;
    input funct7;
    output [2:0] ALUctrl;

    parameter ADD  = 3'd0;
    parameter SUB  = 3'd1;
    parameter AND  = 3'd2;
    parameter XOR  = 3'd3;
    parameter SLT  = 3'd4;
    parameter SRA  = 3'd5;
    parameter SLL  = 3'd6;
    case (ALUop)
        00: ALUctrl = ADD;
        01: ALUctrl = SUB;
        10: case (funct3)
            3'b000: ALUctrl = (funct7) ? SUB:ADD;
            3'b001: ALUctrl = SLL;
            3'b010: ALUctrl = SLT;
            3'b100: ALUctrl = XOR;
            3'b101: ALUctrl = SRA;
            3'b111: ALUctrl = AND;
            default: 
            endcase
        default: 
    endcase

endmodule
module ALU(in1, in2, ALUctrl, result, zero);
    parameter DATA_W = 32;
    input  [DATA_W-1:0] in1; 
    input  [DATA_W-1:0] in2;
    input  [2:0]        ALUctrl;
    output [DATA_W-1:0] result;
    output              zero;

    wire signed [DATA_W-1: 0] signed_in1, signed_in2;
    reg [DATA_W-1 : 0] result_reg;
    reg zero_reg;

    parameter ADD  = 3'd0;
    parameter SUB  = 3'd1;
    parameter AND  = 3'd2;
    parameter XOR  = 3'd3;
    parameter SLT  = 3'd4;
    parameter SRA  = 3'd5;
    parameter SLL  = 3'd6;

    assign signed_in1 = in1;
    assign signed_in2 = in2;
    assign result = result_reg;

    always @(*) begin
        case (ALUctrl)
            ADD: begin
                result_reg[31:0] = in1 + in2;
                //overflow
                if (in1[31] == in2[31]) begin
                    if (~in1[31] && result_reg[31]) result_reg[31:0] = {1'b0,{31{1'b1}}};
                    else if (in1[31] && ~result_reg[31]) result_reg[31:0] = {1'b1,{31{1'b0}}};
                end
            end
            SUB: begin
                result_reg[31:0] = signed_in1 - signed_in2;
                //overflow
                if (in1[31] != in2[31]) begin
                    if (~in1[31] && in2[31] && result_reg[31]) result_reg[31:0] = {1'b0,{31{1'b1}}};
                    else if (in1[31] && ~in2[31] && ~result_reg[31]) result_reg[31:0] = {1'b1,{31{1'b0}}};
                end
                //zero for beq
                if(result_reg == 32'b0) zero = 1'b1;
            end
            AND: result_reg[31:0] = in1 & in2;
            XOR: result_reg[31:0] = in1 ^ in2;
            SLT: result_reg[31:0] = (signed_in1 < signed_in2)? 1:0;
            SRA: result_reg[31:0] = {{32{in1[31]}},in1} >> in2;
            SLL: result_reg[31:0] = in1 << in2;
        endcase
    end
endmodule


module MULDIV_unit(
    result, o_done, i_clk, i_valid, i_A, i_B, ALUctrl
    );
    // Todo: HW2
    parameter DATA_W = 32;
    input  i_clk, i_valid;
    input  [DATA_W - 1 : 0]     i_A;    // input operand A
    input  [DATA_W - 1 : 0]     i_B;    // input operand B
    input  [         2 : 0]     ALUctrl; // instruction
    output [DATA_W - 1 : 0]     result; // output value
    output                      o_done; // output valid signal
//Parameter
    //ALUctrl
    parameter MUL  = 3'd7;
    parameter DIV  = 3'd8;
    //state
    parameter S_IDLE           = 2'd0;
    parameter S_MULTI_CYCLE_OP = 2'd1;

    reg  [         1: 0] state, state_nxt;
    reg  [  DATA_W-1: 0] operand_a, operand_a_nxt;
    reg  [  DATA_W-1: 0] operand_b, operand_b_nxt;
    reg  [         2: 0] inst, inst_nxt;
    reg  [         5: 0] counter;
    reg  [  2*DATA_W: 0] result_reg;
    reg  o_done_reg;

    assign result[DATA_W-1:0] = result_reg[DATA_W-1:0];
    assign o_done = o_done_reg;

    always @(*) begin
        if (i_valid) begin
            operand_a_nxt = i_A;
            operand_b_nxt = i_B;
            inst_nxt      = ALUctrl;
        end
        else begin
            operand_a_nxt = operand_a;
            operand_b_nxt = operand_b;
            inst_nxt      = inst;
        end
    end
    //FSM
    always @(*) begin
        case(state)
            S_IDLE           : begin
                if(!i_valid) state_nxt <= S_IDLE;
                else state_nxt <= (ALUctrl < 3'd7) ? S_ONE_CYCLE_OP : S_MULTI_CYCLE_OP;
            end
            S_MULTI_CYCLE_OP : state_nxt <= (counter == 32) ? S_IDLE : state;
            default : state_nxt <= state;
        endcase
    end
    //counter
    always @(negedge i_clk or negedge i_valid) begin
        if(~i_valid) counter <= counter + 1;
        else counter <= 0;
    end
    //MUL
    always @(negedge i_clk) begin
        if (i_valid) result_reg = i_A;
        if (state == S_MULTI_CYCLE_OP) begin
            case (ALUctrl)
                MUL:begin
                    if(counter >= 0)begin
                        result_reg = (result_reg[0]) ? {result_reg[64:32] + operand_b, result_reg[31:0]}: result_reg;
                        result_reg = result_reg >> 1;
                    end
                    else result_reg = 0;
                end
                DIV: begin
                    if(counter >= 0)begin
                        result_reg = result_reg << 1;
                        result_reg = (result_reg[63:32] >= operand_b) ? {result_reg[64:32] - operand_b, result_reg[31:0]}+1'b1: result_reg;
                    end
                    else result_reg = 0;
                end
                default: result_reg = 0;
            endcase
        end
    end
    //o_done signal
    always @(posedge i_clk) begin
        case (state)
            S_MULTI_CYCLE_OP : o_done_reg <= (counter == 32)? 1: 0;
            default: o_done_reg <= 0;
        endcase
    end
    //Sequential always block
    always @(posedge i_clk or negedge i_rst_n) begin
        if (!i_rst_n) begin
            state       <= S_IDLE;
            operand_a   <= 0;
            operand_b   <= 0;
            inst        <= 0;
        end
        else begin
            state       <= state_nxt;
            operand_a   <= operand_a_nxt;
            operand_b   <= operand_b_nxt;
            inst        <= inst_nxt;
        end
    end
endmodule

module Cache#(
        parameter BIT_W = 32,
        parameter ADDR_W = 32
    )(
        input i_clk,
        input i_rst_n,
        // processor interface
            input i_proc_cen,
            input i_proc_wen,
            input [ADDR_W-1:0] i_proc_addr,
            input [BIT_W-1:0]  i_proc_wdata,
            output [BIT_W-1:0] o_proc_rdata,
            output o_proc_stall,
            input i_proc_finish,
            output o_cache_finish,
        // memory interface
            output o_mem_cen,
            output o_mem_wen,
            output [ADDR_W-1:0] o_mem_addr,
            output [BIT_W*4-1:0]  o_mem_wdata,
            input [BIT_W*4-1:0] i_mem_rdata,
            input i_mem_stall,
            output o_cache_available
    );

    assign o_cache_available = 0; // change this value to 1 if the cache is implemented

    //------------------------------------------//
    //          default connection              //
    assign o_mem_cen = i_proc_cen;              //
    assign o_mem_wen = i_proc_wen;              //
    assign o_mem_addr = i_proc_addr;            //
    assign o_mem_wdata = i_proc_wdata;          //
    assign o_proc_rdata = i_mem_rdata[0+:BIT_W];//
    assign o_proc_stall = i_mem_stall;          //
    //------------------------------------------//

    // Todo: BONUS

endmodule
